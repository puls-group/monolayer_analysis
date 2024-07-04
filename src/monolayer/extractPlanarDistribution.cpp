//#include "Utils.h"
#include "PlanarGeometry.h"
#include "OptionHandler.hpp"
//#include "gmxUtils.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/math/vectypes.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iomanip>
#include <omp.h>

//#define PRECISE_LIFETIMES 0
#define PRINT_LIFETIMES true

// Main
int main(int argc, char* argv[]) {
	using std::vector;
	using std::string;
	using std::cerr;
	using std::flush;
	using std::endl;
	using namespace gmx;
	// Generate logfile
	std::ofstream logfile;
	logfile.open("logfile.txt", std::ofstream::trunc);
	logfile.close();
	logfile.open("logfile.txt", std::ofstream::app);

	std::time_t starttime =
		std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	logfile << "Start date:\n" << std::ctime(&starttime) << std::endl;

	// Parameter potentially modified by parsing
	string input_xtc_1 = "./input1.xtc";
	string input_xtc_2 = "./input2.xtc";
	string output_prefix = "./2d_distr_";
	unsigned long num_mols_1 = 0;
	unsigned long num_mols_2 = 0;
	unsigned long offset_mols_1 = 0;
	unsigned long offset_mols_2 = 0;
	unsigned long distribution_buckets = 201;
	unsigned max_threads = omp_get_max_threads();
	double maxdz = -1;
	bool skipsame = false;

	op::OptionHandler OH(
R"HELP(Tool to calculate 2D relative density distributions of sets of two different populations of vertices/molecules in a monolayer configuration. If not provided with a monolayer configuration, the distribution of the projection onto the `xy`-plain will be calculated. 
It takes two xtc input files/trajectories -i1 and -i2 with the positions of, e.g., cations' and anion' centers of masses respectively but any other separation into two populations is possible. Both must have the same number of frames or the analysis will fail. The tools does not account for a change in the number of molecules per frame within a trajectory. It, however, does not rely on equal numbers of molecules in both trajectories. 

It then proceeds to iterate through the frames in the trajectory, and calculate the distribution of the second population relative to the first population in a grid of configurable size and resolution and it does the same for the distribution of other members of the first population relative to members of the first population (pairs of identical indices will be skipped if  `-s 1` is set). 
The grid is chosen to bin the coefficients of the PBC basis vectors in the `xy`-plain after linear decomposition of relative positions into the basis representation clipped to a `[-0.5, +0.5]` interval. Absolute distributions can be obtained by multiplying with the respective basis vector again.
It furthermore averages its resulting statistics across all frames in the provided trajectories.
In principle, it would also be possible to provide individual atom positions to the code instead of COM for full molecules.
This would, however smear out the resulting distributions especially in `solid-like` systems.

Eventually, this program outputs a set of statistics files into the path denoted by the option `-o`: 
- 1 File `basis.dat` containing the basis vectors in the `xy`-plain used for the basis decomposition. These are extracted from the PBC data associated with the first input trajectory `-i1`
- 1 File `xy_distribution.dat` containing statistics on the distribution of the coeffcients in the basis representation of relative distances associated with the basis vectors in `basis.dat`. 

Please be aware that i1 and i2 may denote the same input trajectory file in xtc format, as for this tool (in contrast to `extractFilmParameters`), the number of molecules of either population as well as the offset of the population in the trajectory can be provided via parameters of the tool**)HELP");

	op::SingleValuePositionalOption<std::string> opt_fni1("input_xtc_1", input_xtc_1, true);
	opt_fni1.description(
		"Path to first input xtc trajectory to be processed.");
	OH.addOption(opt_fni1);
	op::SingleValuePositionalOption<std::string> opt_fni2("input_xtc_2", input_xtc_2, true);
	opt_fni2.description(
		"Path to second input xtc trajectory to be processed.");
	OH.addOption(opt_fni2);

	op::SingleValuePositionalOption<std::string> opt_fno("output_prefix", output_prefix);
	opt_fno.description(
		"Path prefix to output data to be produced.");
	OH.addOption(opt_fno);

	op::SingleValuePositionalOption<unsigned long> opt_nmol1("mols_count_1", num_mols_1, true);
	opt_nmol1.description("Number of molecules of first kind");
	OH.addOption(opt_nmol1);

	op::SingleValuePositionalOption<unsigned long> opt_nmol2("mols_count_2", num_mols_2, true);
	opt_nmol2.description("Number of molecules of second kind");
	OH.addOption(opt_nmol2);

	op::SingleValuePositionalOption<unsigned long> opt_off1("mols_offset_1", offset_mols_1, true);
	opt_off1.description("Offset of molecules of first kind");
	OH.addOption(opt_off1);

	op::SingleValuePositionalOption<unsigned long> opt_off2("mols_offset_2", offset_mols_2, true);
	opt_off2.description("Offset of molecules of second kind");
	OH.addOption(opt_off2);

	op::SingleValueOption<unsigned long> opt_buckets("b", distribution_buckets);
	opt_buckets.description("Number of buckets used for the distribution in any dimension");
	OH.addOption(opt_buckets);

	op::SingleValueOption<double> opt_dz("d", maxdz);
	opt_dz.description("Maximum z-direction between pairs to be considered for distribution. Set <0 for no limit.");
	OH.addOption(opt_dz);

	op::SingleValueOption<bool> opt_skip("s", skipsame);
	opt_skip.description("Flag to skip same index in both kinds.");
	OH.addOption(opt_skip);

	op::pRes res = OH.procOptions(argc, argv);
	if (res != op::ok)
		return res;
	if (maxdz < 0) {
		maxdz = std::numeric_limits<double>::max();
	}

	unsigned distribution_buckets_per_dir = distribution_buckets / 2;

	logfile << flush;

	struct t_fileio *input1, *input2;

	cerr << "Opening first xtc data file " << input_xtc_1 << endl;
	input1 = open_xtc(input_xtc_1.c_str(), "r");
	cerr << "Opening seconf xtc data file " << input_xtc_2 << endl;
	input2 = open_xtc(input_xtc_2.c_str(), "r");


	int natoms1;
	int64_t step1;
	real time1;
	matrix box1;
	rvec* x1 = nullptr;
	real prec1;
	gmx_bool bOK1;
	int cframe = 0;
	v2d base1, base2;
	std::pair<v2d, v2d> pbc;
	
	int natoms2;
	int64_t step2;
	real time2;
	matrix box2;
	rvec* x2 = nullptr;
	real prec2;
	gmx_bool bOK2;

	int exit_code = 0;

	vector<vector<vector<long>>> local_distribution(
		max_threads, vector<vector<long>>(
			2 * distribution_buckets_per_dir + 1 /*dim1*/,
			vector<long>(2 * distribution_buckets_per_dir + 1 /*dim2*/, 0)));

	vector<vector<long>> distribution(
		distribution_buckets /*dim1*/,
		vector<long>(distribution_buckets /*dim2*/, 0));

	std::ofstream file;
	long long total_hits = 0;
	vector<long long> l_total_hits(max_threads, 0);

	printf("\nInitialization successful. Start processing...\n\n");


	for (; (cframe == 0 && read_first_xtc(input1, &natoms1, &step1, &time1, box1, &x1, &prec1, &bOK1) && read_first_xtc(input2, &natoms2, &step2, &time2, box2, &x2, &prec2, &bOK2)) ||
		(read_next_xtc(input1, natoms1, &step1, &time1, box1, x1, &prec1, &bOK1) > 0 && read_next_xtc(input2, natoms2, &step2, &time2, box2, x2, &prec2, &bOK2) > 0);) {
		if (!bOK1 || ! bOK2) {
			cerr << "Frame #" << cframe << " in either i1 or i2 was corrupted" << endl;
			exit_code = 1;
			goto _finalize;
		}

		if (cframe == 0) {
			base1 = v2d(box1[0][0], box1[0][1]);
			base2 = v2d(box1[1][0], box1[1][1]);
			cerr << "Box: " << base1 << " x " << base2 << endl;

			pbc.first = base1;
			pbc.second = base2;
		}
#pragma omp parallel for
		for (unsigned f = 0; f < num_mols_1; f++) {
			int thread_id = omp_get_thread_num();
			//cerr << "r index " << offset_mols_1 + f << endl;
			rvec &rp = x1[offset_mols_1 + f];
			double rz = rp[2];
			v2d rpos(rp[0], rp[1]);
			for (unsigned s = 0; s < num_mols_2; s++) {
				if (skipsame && f == s) {
					continue;
				}
				//cerr << "t index " << offset_mols_2 + s << endl;
				rvec &tp = x2[offset_mols_2 + s];
				double tz = tp[2];
				if (abs(tz - rz) >= maxdz) {
					continue;
				}
				v2d tpos(tp[0], tp[1]);

				//cerr << "rpos: " << rpos << " tpos" << tpos << endl;

				v2d dist = tpos - rpos;
				v2d coeffs = represent(restrict(dist, pbc), pbc);
				//cerr << coeffs << endl;

				double xpos = coeffs.real();
				xpos = std::max(xpos, -0.5);
				xpos = std::min(xpos, 0.5);
				double ypos = coeffs.imag();
				ypos = std::max(ypos, -0.5);
				ypos = std::min(ypos, 0.5);

				unsigned fbuck = round(2.0*xpos*distribution_buckets_per_dir) + distribution_buckets_per_dir;
				unsigned sbuck = round(2.0*ypos*distribution_buckets_per_dir) + distribution_buckets_per_dir;
				//cerr << coeffs << "->" << fbuck << "|" << sbuck << endl;
				local_distribution[thread_id][fbuck][sbuck] ++;
			}
		}

		cframe++;
		if (cframe % 100 == 0) {
			cerr << "\rProcessed frame #" << cframe; // << endl;
		}
	}
	cerr << endl << "Finished processing." << endl;
	cerr << "Aggregating data." << endl;

#pragma omp parallel for
	for (unsigned f = 0; f < 2 * distribution_buckets_per_dir + 1; f++) {
		int thread_id = omp_get_thread_num();
		for (unsigned s = 0; s < 2 * distribution_buckets_per_dir + 1; s++) {
			for (unsigned t = 0; t < max_threads; t++) {
				l_total_hits[thread_id] += local_distribution[t][f][s];
				distribution[f][s] += local_distribution[t][f][s];
			}
		}
	}

	for (unsigned t = 0; t < max_threads; t++) {
		total_hits += l_total_hits[t];
	}

	file.open(output_prefix + "basis.dat",
		std::ofstream::trunc);

	file << "#Each line holds one basis vector" << endl;
	file << "#base x [nm]\t#base y [nm]" << endl;
	file << base1.real() << "\t" << base1.imag() << endl;
	file << base2.real() << "\t" << base2.imag() << endl;
	file.close();

	file.open(output_prefix + "xy_distribution.dat",
		std::ofstream::trunc);

	file << "# The file should be read as a 2 dimensional matrix" << endl;
	file << "# The first line has one dummy entry and the multiples of the y basis associated with the bins in that column " << endl;
	file << "# Each following line has the multiple of the x basis associated with the bins in that row and the mean number of COM in a frame in that bin relative to the total number of COM detected." << endl;
	file << std::setprecision(10) << std::fixed;

	file << 0.0;
	for (unsigned s = 0; s < 2 * distribution_buckets_per_dir + 1; s++) {
		file << "\t" << 0.5*((double(s) - double(distribution_buckets_per_dir))) / distribution_buckets_per_dir;
	}
	file << endl;

	for (unsigned f = 0; f < 2 * distribution_buckets_per_dir + 1; f++) {
		file << 0.5*((double(f) - double(distribution_buckets_per_dir))) / distribution_buckets_per_dir;
		for (unsigned s = 0; s < 2 * distribution_buckets_per_dir + 1; s++) {
			file << "\t" << (double)distribution[f][s] / (double)total_hits;
		}
		file << endl;
	}

	file.close();
_finalize:
	close_xtc(input1);
	close_xtc(input2);
	free(x1);
	free(x2);

	printf("\nStats for xy distribution have been evaluated.\n\n");
	printf("##########################################################\n\n\n");
	std::time_t endtime =
		std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	logfile << "\nEnd date:\n" << std::ctime(&endtime) << std::endl;
	logfile.close();
	return exit_code;
}
