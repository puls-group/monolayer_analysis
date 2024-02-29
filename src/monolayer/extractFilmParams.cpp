#include "OptionHandler.hpp"
#include "GraphConstructor.h"
#include "gromacs/fileio/xtcio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
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

void print(gmx::RVec &pos)
{
	std::cerr << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")"
			  << std::endl;
}

// Main
int main(int argc, char *argv[])
{
	using namespace gmx;
	using namespace std;
	// Parameter potentially modified by parsing
	std::string input_xtc_1 = "./tmp/in1.xtc";
	std::string input_xtc_2 = "./tmp/in2.xtc";
	std::string output_prefix = "/tmp/film_";

	double cutoff_distance = 6;
	unsigned min_hole_size = 6;
	const PbcType pbc_type = PbcType::Xyz;

	bool dump_graphs = false;
	size_t dump_frame_index = 0;

	op::OptionHandler OH(
		R"HELP(Tool to calculate the area and cicumference of monolayer films on sapphire.)HELP");

	op::SingleValueOption<std::string> opt_fni1("i1", input_xtc_1);
	opt_fni1.description("Path to input xtc trajectory to be processed.");
	OH.addOption(opt_fni1);

	op::SingleValueOption<std::string> opt_fni2("i2", input_xtc_2);
	opt_fni2.description("Path to input xtc trajectory to be processed.");
	OH.addOption(opt_fni2);

	op::SingleValueOption<std::string> opt_fno1("o", output_prefix);
	opt_fno1.description("Path to output file of circumference and area data.");
	OH.addOption(opt_fno1);

	op::SingleValueOption<double> opt_cutoff("c", cutoff_distance);
	opt_cutoff.description(
		"Cutoff distance up to which the two particles should be considered connected in graph.");
	OH.addOption(opt_cutoff);

	op::SingleValueOption<unsigned> opt_hole_lower("hl", min_hole_size);
	opt_hole_lower.description(
		"Lowest number of corners a face must have to be considered a hole in the film.");
	OH.addOption(opt_hole_lower);

	op::SingleValueOption<bool> opt_dump_graphs("d", dump_graphs);
	opt_dump_graphs.description(
		"Dump the first frame graph representation and exit.");
	OH.addOption(opt_dump_graphs);
	op::SingleValueOption<size_t> opt_dump_frame_index("df", dump_frame_index);
	opt_dump_frame_index.description(
		"Index of the frame to dump the graph of.");
	OH.addOption(opt_dump_frame_index);

	op::pRes res = OH.procOptions(argc, argv);
	if (res != op::ok)
		return res;
	struct t_fileio *input1, *input2;

	// double cutoff_distance_sq = cutoff_distance * cutoff_distance;
	input1 = open_xtc(input_xtc_1.c_str(), "r");
	input2 = open_xtc(input_xtc_2.c_str(), "r");

	ofstream output;
	ofstream stats_output;
	ofstream holes_output;
	std::ofstream graph_out;
	std::ofstream points_out;

	int exit_code = 0;

	// unsigned int max_threads = omp_get_max_threads();
	int natoms1, natoms2;
	int64_t step1, step2;
	::real time1, time2;
	matrix box1, box2;
	rvec *x1, *x2;
	::real prec1, prec2;
	gmx_bool bOK;
	int cframe = 0;
	t_pbc pbc_info;
	vector<v2d> graph_nodes;
	pair<v2d, v2d> graph_pbc;
	Graph2D *graph;

	vector<size_t> hole_size_statistics;
	vector<size_t> hole_count_statistics;
	double hole_count_accumulated = 0;

	size_t total_node_count;
	size_t total_stats_count;
	size_t iteration_limit;
	std::pair<std::vector<size_t>, size_t> stats;
	std::vector<std::pair<size_t, std::pair<double, double>>> hole_stats;

	if (!read_first_xtc(input1, &natoms1, &step1, &time1, box1, &x1, &prec1,
						&bOK))
	{
		cerr << "Failed to read from input 1" << endl;
		exit_code = 1;
		goto _finalize;
	}
	if (!bOK)
	{
		cerr << "Frame was corrupted in input 1" << endl;
		exit_code = 1;
		goto _finalize;
	}
	if (!read_first_xtc(input2, &natoms2, &step2, &time2, box2, &x2, &prec2,
						&bOK))
	{
		cerr << "Failed to read from input 2" << endl;
		exit_code = 1;
		goto _finalize;
	}
	if (!bOK)
	{
		cerr << "Frame was corrupted in input 2" << endl;
		exit_code = 1;
		goto _finalize;
	}

	// Debug dump of the PBC conditions
	cerr << box1[0][0] << "\t" << box1[0][1] << "\t" << box1[0][2] << endl;
	cerr << box1[1][0] << "\t" << box1[1][1] << "\t" << box1[1][2] << endl;
	cerr << box1[2][0] << "\t" << box1[2][1] << "\t" << box1[2][2] << endl;

	set_pbc(&pbc_info, pbc_type, box1);

	total_node_count = natoms1 + natoms2;
	total_stats_count = total_node_count * 2 + 1;
	// Allocate all required memory
	graph_nodes.resize(total_node_count);
	hole_size_statistics.resize(total_stats_count);
	hole_count_statistics.resize(total_stats_count);

	graph_pbc = {{box1[0][0], box1[0][1]}, {box1[1][0], box1[1][1]}};
	// Iterate over first kind of ions
#pragma omp parallel for
	for (unsigned f = 0; f < unsigned(natoms1); f++)
	{
		graph_nodes[f].real(x1[f][0]);
		graph_nodes[f].imag(x1[f][1]);
	}
	// Add other kind of ions
#pragma omp parallel for
	for (unsigned s = 0; s < unsigned(natoms2); s++)
	{
		graph_nodes[natoms1 + s].real(x2[s][0]);
		graph_nodes[natoms1 + s].imag(x2[s][1]);
	}

	graph = constructPlanarGraph(graph_nodes, graph_pbc, cutoff_distance);
	stats = graph->findFaceStatistics(min_hole_size);
	hole_stats = graph->findHoleStatistics(min_hole_size);

	if (dump_graphs && dump_frame_index == cframe)
	{
		graph_out.open(output_prefix + "graph_dump.svg", ofstream::trunc);
		graph->dump(graph_out, min_hole_size);
		graph_out.close();

		points_out.open(output_prefix + "points.txt", ofstream::trunc);
		graph->dumpNodes(points_out);
		points_out.close();
		goto _finalize;
	}

	output.open(output_prefix + "outline_area.dat", ofstream::trunc);
	output << "#Frame Nr [1]\t#circumference[nm]\t#area[nm^2]\t#number_of_holes[1]" << endl;
	output << cframe << "\t" << graph->findOutline(min_hole_size) << "\t" << graph->findArea(min_hole_size) << "\t" << stats.second << endl;

	holes_output.open(output_prefix + "holes_statistics.dat", ofstream::trunc);
	holes_output << "#List of all holes detected in the film trajectory\n";
	holes_output << "# Number of vertices\t#circumference[nm]\t#area[nm^2]" << endl;

	holes_output << "#" << cframe << "\n";

	for (const std::pair<size_t, std::pair<double, double>> &entry : hole_stats)
	{
		holes_output << entry.first << "\t" << entry.second.first << "\t" << entry.second.second << "\n";
	}

	holes_output.flush();

	hole_count_statistics[std::min(hole_count_statistics.size() - 1, stats.second)]++;
	hole_count_accumulated += stats.second;

	iteration_limit = std::min(total_stats_count, stats.first.size());
	for (size_t i = 0; i < iteration_limit; i++)
	{
		hole_size_statistics[i] += stats.first[i];
	}

	delete graph;

	cframe++;

	for (; read_next_xtc(input1, natoms1, &step1, &time1, box1, x1, &prec1,
						 &bOK) > 0;)
	{
		if (!bOK)
		{
			cerr << "Frame was corrupted in input 1" << endl;
			exit_code = 1;
			goto _finalize;
		}

		read_next_xtc(input2, natoms2, &step2, &time2, box2, x2, &prec2, &bOK);
		if (!bOK)
		{
			cerr << "Frame was corrupted in input 2" << endl;
			exit_code = 1;
			goto _finalize;
		}

		// Iterate over first kind of ions
#pragma omp parallel for
		for (unsigned f = 0; f < unsigned(natoms1); f++)
		{
			graph_nodes[f].real(x1[f][0]);
			graph_nodes[f].imag(x1[f][1]);
		}
		// Add other kind of ions
#pragma omp parallel for
		for (unsigned s = 0; s < unsigned(natoms2); s++)
		{
			graph_nodes[natoms1 + s].real(x2[s][0]);
			graph_nodes[natoms1 + s].imag(x2[s][1]);
		}

		graph = constructPlanarGraph(graph_nodes, graph_pbc, cutoff_distance);

		stats = graph->findFaceStatistics(min_hole_size);
		hole_stats = graph->findHoleStatistics(min_hole_size);

		if (dump_graphs && dump_frame_index == cframe)
		{
			graph_out.open(output_prefix + "graph_dump.svg", ofstream::trunc);
			graph->dump(graph_out, min_hole_size);
			graph_out.close();

			points_out.open(output_prefix + "points.txt", ofstream::trunc);
			graph->dumpNodes(points_out);
			points_out.close();
			goto _finalize;
		}

		output << cframe << "\t" << graph->findOutline(min_hole_size) << "\t" << graph->findArea(min_hole_size) << "\t" << stats.second << endl;

		holes_output << "#" << cframe << "\n";

		for (const std::pair<size_t, std::pair<double, double>> &entry : hole_stats)
		{
			holes_output << entry.first << "\t" << entry.second.first << "\t" << entry.second.second << "\n";
		}

		holes_output.flush();

		hole_count_statistics[std::min(hole_count_statistics.size() - 1, stats.second)]++;
		hole_count_accumulated += stats.second;

		size_t iteration_limit = std::min(total_stats_count, stats.first.size());
		for (size_t i = 0; i < iteration_limit; i++)
		{
			hole_size_statistics[i] += stats.first[i];
		}
		delete graph;

		cframe++;
		if (cframe % 100 == 0)
		{
			cerr << "\rAnalyzed frame #" << cframe; // << endl;
		}
	}
	cerr << endl;
	cerr << "Saving hole statistics" << endl;

	stats_output.open(output_prefix + "face_statistics.dat", ofstream::trunc);

	stats_output << "# This file has the statistics of mean holes of a certain size per frame" << endl;
	stats_output << "# The film on average had " << std::scientific << std::setprecision(4) << (double(hole_count_accumulated) / double(cframe)) << " holes per frame" << endl;
	stats_output << "# Up to a maximum of " << (total_stats_count - 1) << " vertices per hole, the distribution was as follows: " << endl;
	stats_output << "# number of vertices\t#number of mean occurrences per frame" << endl;
	stats_output << std::scientific << std::setprecision(8);
	for (size_t i = 0; i < iteration_limit; i++)
	{
		if (hole_size_statistics[i] > 0)
		{
			stats_output << i << "\t" << (double(hole_size_statistics[i]) / double(cframe)) << endl;
		}
	}
	stats_output.close();
	cerr << "Analysis done. Cleaning up..." << endl;
_finalize:
	close_xtc(input1);
	close_xtc(input2);
	free(x1);
	free(x2);
	output.close();
	holes_output.close();

	return exit_code;
}
