#include "GraphConstructor.h"
#include <iostream>

using namespace std;

/*
int main(int argc, char**argv) {
	Graph2D graph(5);

	vector<v2d> positions = { {0.0,0.0}, {0.0, 2.0}, {1, 0.0}, {2,0.0}, {1.5,4} };
	for (int i = 0; i < 5; i++) {
		graph.getNode(i)->setPosition(positions[i]);
	}

	graph.connect(0, 1, positions[1] - positions[0]);
	cerr << abs(positions[1] - positions[0]) << endl;
	graph.connect(0, 2, positions[2] - positions[0]);
	cerr << abs(positions[2] - positions[0]) << endl;
	graph.connect(0, 4, positions[4] - positions[0]);
	cerr << abs(positions[4] - positions[0]) << endl;
	graph.connect(2, 3, positions[3] - positions[2]);
	cerr << abs(positions[3] - positions[2]) << endl;
	graph.connect(3, 4, positions[4] - positions[3]);
	cerr << abs(positions[4] - positions[3]) << endl;
	graph.connect(2, 4, positions[4] - positions[2]);
	cerr << abs(positions[4] - positions[2]) << endl;

	cerr << graph.findOutline(5) << endl;
	graph.dump(cout);
	cout << endl;

	graph.reinit(5);
	for (int i = 0; i < 5; i++) {
		graph.getNode(i)->setPosition(positions[i]);
	}

	graph.connect(0, 1, positions[1] - positions[0]);
	cerr << abs(positions[1] - positions[0]) << endl;
	graph.connect(0, 2, positions[2] - positions[0]);
	cerr << abs(positions[2] - positions[0]) << endl;
	graph.connect(2, 3, positions[3] - positions[2]);
	cerr << abs(positions[3] - positions[2]) << endl;
	graph.connect(3, 4, positions[4] - positions[3]);
	cerr << abs(positions[4] - positions[3]) << endl;
	graph.connect(4, 1, positions[1] - positions[4]);
	cerr << abs(positions[1] - positions[4]) << endl;

	cerr << graph.findOutline(15) << endl;

	graph.dump(cout);
	cout << endl;
	return 0;
}*/

int main(int argc, char**argv) {

	pair<v2d, v2d> pbc = { { 0,3 },{ 3,0 } };

	v2d testv(4, 0);
	cerr << "Check basis decomposition of vector" << endl;
	v2d coeff = represent(testv, pbc);
	cerr << testv << "-->" << pbc.first << "++" << pbc.second << "::" << coeff << endl;
	coeff = restrict(testv, pbc);
	cerr << testv << "-->" << pbc.first << "++" << pbc.second << "::" << coeff << endl;

	cerr << "Check pbc graph construction" << endl;
	vector<v2d> points = { { 0.5,0.5 },{ 0.5,1.5 },{ 0.5,2.5 },{ 1.5,2 },{ 2.5,2.5 },{ 2,1 } };

	Graph2D *res = constructPlanarGraph(points, pbc, sqrt(2.55));
	cerr << "done with graph creation" << endl;
	cerr << "Check graph analysis routines" << endl;
	res->findFaces();
	res->dump(cout);
	cout << endl;
	cerr << "Done dumping" << endl;

	cerr << "Outline: " << res->findOutline(1000) << endl;
	cerr << "Area: " << res->findArea(1000) << endl;
	//res->dump(cout);
	//cout << endl;
	delete res;


	pbc = { { 0,100},{ 100,0 } };

	points = { { 0.5,0.5 },{ 0.5,1.5 },{ 0.5,2.5 },{ 1.5,2 },{ 2.5,2.5 },{ 2,1 } };

	cerr << "Check finite graph analysis" << endl;
	res = constructPlanarGraph(points, pbc, sqrt(2.55));
	cerr << "done with graph creation" << endl;
	res->findFaces();
	res->dump(cout);
	cout << endl;
	cerr << "Done dumping" << endl;

	cerr << "Outline: " << res->findOutline(1000) << endl;
	cerr << "Area: " << res->findArea(1000) << endl;
	//res->dump(cout);
	//cout << endl;
	delete res;



	return 0;
}
