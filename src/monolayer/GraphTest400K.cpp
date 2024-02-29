#include "GraphConstructor.h"
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv)
{

	pair<v2d, v2d> pbc = {{46.5462, 0.}, {-7.7577, 13.4364}};

	vector<v2d> points;

	std::ifstream points_in;
	std::ofstream graph_out,graph_out2;

	points_in.open("./points.txt", std::ios::in);
	graph_out.open("./graph_400K.svg", std::ios::out);
	graph_out2.open("./graph_intersection.svg", std::ios::out);

	size_t num_points;
	points_in >> num_points;

	for(size_t p = 0; p < num_points; p++){
		double real_, imag_;

		points_in >> real_ >> imag_;
		points.push_back({real_, imag_});
	}

	size_t min_hole_size = 6;

	Graph2D *res;
	res = constructPlanarGraph(points, pbc, sqrt(2));
	cerr << "done with graph creation" << endl;
	cerr << "Check graph analysis routines" << endl;
	res->findFaces();
	res->dump(graph_out, min_hole_size);
	cout << endl;
	cerr << "Done dumping" << endl;
	delete res;

	points = {{0,1},{1,0},{0,0},{1,1}};
	res = constructPlanarGraph(points, pbc, sqrt(2));
	res->findFaces();
	res->dump(graph_out2, min_hole_size);
	cerr << "Done dumping intersection graph" << endl;
	delete res;


	return 0;
}
