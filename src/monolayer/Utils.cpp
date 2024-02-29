#include "Utils.h"
//#include "gmxcpp/Trajectory.h"
//#include "gmxcpp/Utils.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

MSD_History::MSD_History(long long max_size) : max_size(max_size) {
    if (max_size >= 0) {
        limited = true;
        history.assign(max_size, {0, 0, 0});
    } else {
        limited = false;
        history.clear();
    }
    curr_start = 0;
    curr_next = curr_start;
    curr_length = 0;
}

int MSD_History::size() {
    if (limited) {
        return curr_length;
    } else {
        return history.size();
    }
}

coordinates MSD_History::getBack(int steps) {
    if (limited) {
        int pos = curr_next - 1 - steps;
        if (pos < 0) {
            pos %= max_size; // fix periodicity
            pos += max_size; // adjust negative space to positive
            pos %= max_size; // fix periodicity again
        }
        return history[pos];
    } else {
        return history[curr_length - 1 - steps];
    }
}

void MSD_History::pushBack(coordinates coord) {
    if (limited) {
        history[curr_next++] = coord;
        curr_next %= max_size;
        if (curr_length < max_size) {
            curr_length++;
        } else {
            curr_start = (curr_start + 1) % max_size;
        }
    } else {
        history.push_back(coord);
    }
}

void MSD_History::clear() {
    if (limited) {
        curr_start = 0;
        curr_next = curr_start;
    } else {
        history.clear();
    }
    curr_length = 0;
}

void printLifetimeData(string& outputpath, string suffix,
                       vector<vector<long long>>& local_lifetime) {

    int binNr = local_lifetime.size();

    for (int bin = 0; bin < binNr; bin++) {
        std::ofstream lt_file;
        lt_file.open(std::string(outputpath) + PathSeparator + "bin_" +
                         std::to_string(bin) + suffix,
                     std::ofstream::trunc);
        lt_file << "# List of lifetimes in this bin" << endl;
        for (long long lt : local_lifetime[bin]) {
            lt_file << lt << endl;
        }
        lt_file.close();
    }
}

void calculate_MSD(vector<double>& msd_pll, vector<double>& msd_perp,
                   vector<long long>& msd_count, MSD_History& position_history,
                   coordinates position) {
    int total_prev = position_history.size();
    for (int dt = 0; dt < total_prev; dt++) {
        coordinates ref = position_history.getBack(dt);
        /*coordinates pll_ref(ref);
        pll_ref[2] = 0.0;
        coordinates perp_ref(ref);
        perp_ref[0] = 0.0;
        perp_ref[1] = 0.0;
        coordinates pll_pos(position);
        pll_pos[2] = 0.0;
        coordinates perp_pos(position);
        perp_pos[0] = 0.0;
        perp_pos[1] = 0.0;
        // Calculate parallel distance squared
        double res = distance2(pll_ref, pll_pos);
        if ((unsigned)dt >= msd_pll.size()) {
            cerr << "dt " << dt << " sz " << msd_pll.size() << endl;
        }
        assert((unsigned)dt < msd_pll.size());
        msd_pll[dt] += res;

        // Calculate perpendicular distance squared
        res = distance2(perp_ref, perp_pos);
        msd_perp[dt] += res;*/
        double dx = ref[0] - position[0];
        double dy = ref[1] - position[1];
        double dz = ref[2] - position[2];
        msd_pll[dt] += dx * dx + dy * dy;
        msd_perp[dt] += dz * dz;

        // Count occurrence of particle in bin
        msd_count[dt]++;
    }

    position_history.pushBack(position);
}

// Version: 3.2
// Hotfix: Fixed survival probability computation

int getAbsoluteBin(std::set<std::pair<float, int>>& binBoundaries,
                   coordinates coord) {
    float z = coord[2];
    return getAbsoluteBin(binBoundaries, z);
}

int getAbsoluteBin(set<pair<float, int>>& binBoundaries, float z) {
    int max = binBoundaries.size();
    auto it = binBoundaries.lower_bound({z, max});

    if (it == binBoundaries.end()) {
        return max;
    } else {
        return it->second;
    }
}

int getBin(std::set<std::pair<float, int>>& absBoundaries,
           std::vector<std::pair<float, float>>& graceBoundaries,
           int previous_bin, coordinates coord) {

    int max = absBoundaries.size();
    float z = coord[2];

    if (previous_bin >= 0 && previous_bin < max) {
        std::pair<float, float> limits = graceBoundaries[previous_bin];
        if (limits.xx <= z && z <= limits.yy) {
            return previous_bin;
        }
    }

    return getAbsoluteBin(absBoundaries, coord);
}

bool isInInterval(float lower, float higher, float val) {
    if (lower > higher) {
        return val <= lower || val >= higher;
    } else {
        return val >= lower && val <= higher;
    }
}

int getPeriodicAbsoluteBin(std::set<std::pair<float, int>>& binBoundaries,
                           coordinates coord,
                           vector<pair<float, float>>& periods) {
    coord = getPeriodCorrectedCoords(coord, periods);
    return getAbsoluteBin(binBoundaries, coord);
}

int getPeriodicBin(set<pair<float, int>>& absBoundaries,
                   vector<pair<float, float>>& graceBoundaries,
                   int previous_bin, coordinates coord,
                   vector<pair<float, float>>& periods) {
    coord = getPeriodCorrectedCoords(coord, periods);
    int max = absBoundaries.size();
    float z = coord[2];

    if (previous_bin >= 0 && previous_bin < max) {
        pair<float, float> limits = graceBoundaries[previous_bin];
        if (isInInterval(limits.xx, limits.yy, z)) {
            return previous_bin;
        }
    }

    return getAbsoluteBin(absBoundaries, coord);
}

coordinates getPeriodCorrectedCoords(coordinates coord,
                                     vector<pair<float, float>>& periods) {
    for (int i = 0; i < 3; i++) {
        pair<float, float> period = periods[i];
        if (period.xx == 0.0f) {
            continue;
        }
        float period_low = period.xx;
        float period_high = period.yy;
        float period_length = period_high - period_low;

        float pos = coord[i];

        float relative = pos - period_low;
        float period_count = floor(relative / period_low);
        coord[i] = pos - period_count * period_length;
    }
    return coord;
}
