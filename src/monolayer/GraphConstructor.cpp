#include "GraphConstructor.h"
#include <iostream>
#include <set>

using namespace std;

Graph2D *constructPlanarGraph(const vector<v2d> &points_, const pair<v2d, v2d> &pbc, double cutoff_dist)
{
	Graph2D *graph = new Graph2D(points_.size());
	using edge_t = pair<pair<int, int>, v2d>;

	vector<edge_t> selected_edges;
	vector<pair<double, edge_t>> shortest_edges;

	vector<v2d> points(points_.size());

	double cutoff_sq = cutoff_dist * cutoff_dist;

	for (unsigned i = 0; i < points.size(); i++)
	{
		points[i] = restrict(points_[i]-pbc.first/2. - pbc.second/2., pbc);
		graph->getNode(i)->setPosition(points[i]);
	}

	// TODO: check length calculation. There seem to be some edges missing.
	vector<v2d> periodic_offsets = {{0, 0}, pbc.first, pbc.second, -pbc.first, -pbc.second, pbc.first + pbc.second, -pbc.first + pbc.second, pbc.first - pbc.second, -pbc.first - pbc.second};
	for (unsigned i = 0; i < points.size(); i++)
	{
		for (unsigned j = 0; j < i; j++)
		{
			v2d dist = points[j] - points[i];
			v2d per_dist = restrict(dist, pbc);
			// cerr << i << "->" << j << ": " << dist << " || " << per_dist << endl;
			double len2 = norm(per_dist);
			if (len2 < cutoff_sq)
			{
				for (const v2d &per : periodic_offsets)
				{
					v2d transvec = per_dist + per;
					len2 = norm(transvec);
					if (len2 < cutoff_sq)
					{
						shortest_edges.push_back({abs(transvec), {{i, j}, transvec}});
					}
				}
			}
			else
			{
				// cerr << "skipped" << endl;
			}
		}
	}
	// cerr << "Total: " << shortest_edges.size() << " edges considered" << endl;

	sort(shortest_edges.begin(), shortest_edges.end(),
		 [](const pair<double, edge_t> &a, const pair<double, edge_t> &b)
		 {
			 return a.first < b.first;
		 });

	/*cerr << "edges sorted:" << endl;
	for (const pair<double, edge_t> &edge : shortest_edges) {
		cerr << edge.first << ": " << edge.second.first.first << "->" << edge.second.first.second << " " << edge.second.second << endl;
	}*/

	for (const pair<double, edge_t> &edge : shortest_edges)
	{
		const edge_t &edat = edge.second;
		int ind1 = edat.first.first;
		int ind2 = edat.first.second;
		bool intersects_previous = false;

		// First option for position of edge in graph
		v2d a = points[ind1];
		v2d b = a + edat.second;

		// Alternate position of edge in graph
		v2d a_ = points[ind2];
		v2d b_ = a_ - edat.second;
		for (edge_t &other : selected_edges)
		{
			int oind1 = other.first.first;
			int oind2 = other.first.second;

			// If endpoints are shared, we do not need to check intersection
			if (ind1 == oind1 || ind1 == oind2 || ind2 == oind1 || ind2 == oind2)
			{
				continue;
			}

			// No shared endpoints, check intersection

			v2d c = points[oind1];
			v2d d = c + other.second;

			// alternate position of edge
			v2d e = points[oind2];
			v2d f = e - other.second;

			for (v2d &offset : periodic_offsets)
			{
				if (intersectsSegment(a, b, c + offset, d + offset) || intersectsSegment(a_, b_, c + offset, d + offset) || intersectsSegment(a, b, e + offset, f + offset) || intersectsSegment(a_, b_, e + offset, f + offset))
				{
					intersects_previous = true;
					goto _quitintersectloop;
				}
			}
		}
	_quitintersectloop:
		if (!intersects_previous)
		{
			selected_edges.push_back(edge.second);
			// cerr << "Longest edge: " << sqrt(edge.first) << " " << endl;
		}
	}

	// cerr << "Checked for intersections, got " << selected_edges.size() << " edges " << endl;

	for (edge_t &other : selected_edges)
	{
		graph->connect(other.first.first, other.first.second, other.second);
	}
	// cerr << "done with graph creation" << endl;
	return graph;
}

/*int main(int argc, char**argv) {

	pair<v2d, v2d> pbc = { { 0,3 },{ 3,0 } };

	vector<v2d> points = { { 0.5,0.5 },{ 0.5,1.5 },{ 0.5,2.5 },{ 1.5,2 },{ 2.5,2.5 },{ 2,1 } };

	Graph2D *res = constructPlanarGraph(points, pbc, sqrt(2.55));
	cerr << "done with graph creation" << endl;
	res->dump(cout);
	cout << endl;
	cerr << "Done dumping" << endl;

	cerr << "Outline: " << res->findOutline() << endl;
	res->dump(cout);
	cout << endl;
	delete res;
	return 0;
}*/