#ifndef __GMXUTILS_H__
#define __GMXUTILS_H__

#include "gromacs/fileio/tpxio.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "TopologyFile.hpp"
#include <map>
#include <set>
#include <string>
#include <vector>
#include <cassert>

bool hasEnding(std::string const &fullString, std::string const &ending) {
	if (fullString.length() >= ending.length()) {
		return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
	}
	else {
		return false;
	}
}

double *read_weights_from_file(std::string file_path, int *tpr_natoms);

/*
* Reads the list of weights from the tpr file and returns an array holding those values.
* The array needs to be deleted after use.
*/
double *read_weights_from_tpr(std::string tpr_path, int *tpr_natoms) {
	assert(tpr_natoms != nullptr);
	struct gmx_mtop_t topology;
	matrix box;
	read_tpx(tpr_path.c_str(), nullptr, box, tpr_natoms, nullptr, nullptr, &topology);
	gmx_mtop_finalize(&topology);

	gmx_mtop_atomloop_all_t iterator = gmx_mtop_atomloop_all_init(&topology);
	t_atom *atom;
	int global_index;

	double *weights = new double[*tpr_natoms];

	while (gmx_mtop_atomloop_all_next(iterator, &global_index, (const t_atom**)&atom)) {
		weights[global_index] = atom->m;
	}
	return weights;
}

double *read_weights_from_wgt(std::string wgt_path, int *tpr_natoms) {
	assert(tpr_natoms != nullptr);

	using namespace topology;

	TopologyWeightFile weightfile(wgt_path, false);
	TopologyHeader *header = weightfile.get_header();
	*tpr_natoms = header->num_entries;

	TopologyWeightEntry entry;
	weightfile.goto_data_start();
	double *res = new double[*tpr_natoms];
	for (unsigned n = 0; n < unsigned(*tpr_natoms); n++) {
		weightfile.read_entry(entry);
		res[n] = entry.weight;
	}

	weightfile.close();

	return res;
}

double *read_weights_from_file(std::string file_path, int *tpr_natoms) {
	if (hasEnding(file_path, ".tpr")) {
		return read_weights_from_tpr(file_path, tpr_natoms);
	}
	else if (hasEnding(file_path, ".wgt")) {
		return read_weights_from_wgt(file_path, tpr_natoms);
	}
	else {
		*tpr_natoms = 0;
		return nullptr;
	}
}

#endif