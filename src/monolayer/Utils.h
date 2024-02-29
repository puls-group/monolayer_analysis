#ifndef __UTILS_H__
#define __UTILS_H__
//#include "gmxcpp/Trajectory.h"
//#include "gmxcpp/Utils.h"
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#define xx first
#define yy second

using std::cerr;
using std::cout;
using std::endl;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

const char PathSeparator =
#if defined _WIN32 || defined __CYGWIN__
    '\\';
#else
    '/';
#endif

struct coordinates {
    float pos[3];

    float& operator[](size_t i) {
        assert(i < 3);
        return pos[i];
    }

    float operator[](size_t i) const {
        assert(i < 3);
        return pos[i];
    }
};

class MSD_History {
  private:
    long long max_size;
    vector<coordinates> history;
    bool limited;
    int curr_start, curr_length, curr_next;

  public:
    MSD_History(long long max_size);
    int size();
    coordinates getBack(int steps);
    void pushBack(coordinates coord);
    void clear();
};

void printLifetimeData(string& outputpath, string suffix,
                       vector<vector<long long>>& lifetime);

void calculate_MSD(vector<double>& msd_pll, vector<double>& msd_perp,
                   vector<long long>& msd_count, MSD_History& position_history,
                   coordinates position);

int getAbsoluteBin(set<pair<float, int>>& binBoundaries, coordinates coord);

int getAbsoluteBin(set<pair<float, int>>& binBoundaries, float z);

int getBin(set<pair<float, int>>& absBoundaries,
           vector<pair<float, float>>& graceBoundaries, int previous_bin,
           coordinates coord);

coordinates getPeriodCorrectedCoords(coordinates coord,
                                     vector<pair<float, float>>& periods);

int getPeriodicAbsoluteBin(set<pair<float, int>>& binBoundaries,
                           coordinates coord,
                           vector<pair<float, float>>& periods);

int getPeriodicBin(set<pair<float, int>>& absBoundaries,
                   vector<pair<float, float>>& graceBoundaries,
                   int previous_bin, coordinates coord,
                   vector<pair<float, float>>& periods);
#endif
