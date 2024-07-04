#include "GraphConstructor.h"
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{

	std::pair<v2d, v2d> pbc = {{46.5462, 0.}, {-7.7577, 13.4364}};

	std::vector<v2d> points;

	std::ifstream points_in;
	std::ofstream graph_out, graph_out2;

	std::string point_input_path = "./tests/points.txt";

	points_in.open(point_input_path, std::ios::in);

	if (points_in.bad())
	{
		std::cerr << "Could not find input file " << point_input_path << " please make sure the test data is present";
		exit(1);
	}

	std::string graph_output_path = "./tests/graph_400K.svg";

	std::cout << "This program reads sample data of a frame of the 400K system of the associated paper and constructs the molecular graph on it." << std::endl;
	std::cout << "The resulting graph is written as an SVG representation to " << graph_output_path << "." << std::endl;
	std::cout << "References for how it should look are provided in ./tests/references/graph_400K.svg." << std::endl;

	graph_out.open(graph_output_path, std::ios::out);

	size_t num_points;
	points_in >> num_points;

	for (size_t p = 0; p < num_points; p++)
	{
		double real_, imag_;

		points_in >> real_ >> imag_;
		points.push_back({real_, imag_});
	}

	size_t min_hole_size = 6;

	Graph2D *res;
	res = constructPlanarGraph(points, pbc, sqrt(2));
	std::cout << "Done with graph creation" << std::endl;
	std::cout << "Checking graph analysis routines" << std::endl;
	res->findFaces();
	res->dump(graph_out, min_hole_size);
	std::cout << std::endl;
	std::cerr << "Done dumping to " << graph_output_path << std::endl;
	delete res;

	return 0;
}
