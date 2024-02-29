#ifndef __GRAPH2D_H__
#define __GRAPH2D_H__
#include "Graph.h"
#include "PlanarGeometry.h"
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <cassert>
// #include <iostream>
#include <complex>
#include <queue>
#include <ostream>
#include <limits>

class Node2D;
class Edge2D;
class Graph2D;

class Edge2D : public Edge
{
protected:
	v2d translation;
	double direction;
	double length;

public:
	Edge2D() : Edge(), translation(0.0, 0.0), direction(0.0), length(0.0) {}

	void setTranslation(double x, double y)
	{
		translation.imag(y);
		translation.real(x);
		direction = std::arg(translation);
		length = std::abs(translation);
	}

	void setTranslation(v2d trans)
	{
		translation = trans;
		direction = std::arg(translation);
		length = std::abs(translation);
	}

	const v2d &getTranslation() const
	{
		return translation;
	}

	double getLength() const
	{
		return length;
	}

	double getDirection() const
	{
		return direction;
	}
};

class Node2D : public Node<Edge2D>
{
protected:
	v2d position;
	// std::vector<std::pair<int, v2d>, int> indexlookup;
public:
	Node2D() : Node<Edge2D>() {}
	virtual ~Node2D() {}

	void setPosition(double x, double y)
	{
		position.imag(y);
		position.real(x);
	}

	void setPosition(v2d pos)
	{
		position = pos;
	}

	virtual void addEdge(Edge2D *edge)
	{
		Node<Edge2D>::addEdge(edge);
		// indexmap[{edge->getTarget(), -edge->getTranslation()}] = this->_adjacent.size() - 1;
	}

	const v2d &getPosition() const
	{
		return position;
	}

	void sortEdges()
	{
		std::sort(_adjacent.begin(), _adjacent.end(),
				  [](const Edge2D *a, const Edge2D *b)
				  {
					  return a->getDirection() > b->getDirection();
				  });

		/*for (unsigned i = 0; i < _adjacent.size(); i++) {
			Edge2D* edge = _adjacent[i];
			//indexmap[edge->getTarget()] = i;
		}*/
	}

	/**
	 * @brief Find the outgoing edge from the current node opposite of the edge with translation verctor trans and pointing towards the point origin.
	 *
	 * origin and trans identify the edge left of which the new edge should be.
	 * I suppose "left of" means next counterclockwise based on how we sort outgoing edges by their argument.
	 *
	 * @param origin Node towards which the edge needs to point
	 * @param trans Opposite translation that the edge needs to have
	 * @return Edge2D*
	 */
	Edge2D *findLeftof(unsigned origin, v2d trans)
	{
		unsigned index = 0;
		for (unsigned i = 0; i < _adjacent.size(); i++)
		{
			Edge2D *e = _adjacent[i];
			if (e->getTarget() == origin && std::norm(trans + e->getTranslation()) < 1e-5)
			{
				index = i;
				break;
			}
		}
		unsigned leftindex = (index + 1) % _adjacent.size();

		if (leftindex >= _adjacent.size() || leftindex < 0)
		{
			return nullptr;
		}

		return _adjacent[leftindex];
	}
};

class Face2D
{
public:
	std::vector<Edge2D *> edges;
	std::vector<int> corners;
	bool isInner;
};

/**
 * @brief Graph class to deal with two dimensional geometric graph in potentially periodic boundary conditions.
 *
 */
class Graph2D : public Graph<Node2D, Edge2D>
{
protected:
	/**
	 * @brief Flags to signify a vertex is a leaf
	 *
	 */
	std::vector<bool> isLeaf;
	/**
	 * @brief Flags to signify a vertex is a branch connecting a leaf to the bulk/inner vertices
	 *
	 */
	std::vector<bool> isBranch;
	/**
	 * @brief Flags to signify a vertex is neither leaf or branch so in the bulk of the connected graph.
	 *
	 */
	std::vector<bool> isInner;

	/**
	 * @brief Flag to mark whether we have performed the face search yet
	 *
	 */
	bool facesfound = false;
	/**
	 * @brief List of found faces in the planar graph
	 *
	 */
	std::vector<Face2D> faces;

	// double circumference;
public:
	Graph2D(int num_nodes) : Graph<Node2D, Edge2D>(num_nodes),
							 isLeaf(num_nodes, false),
							 isBranch(num_nodes, false),
							 isInner(num_nodes, false),
							 facesfound(false) //,
											   // circumference(0.0)
	{
	}

	virtual ~Graph2D() {}

	virtual void reinit(int num_nodes)
	{
		Graph<Node2D, Edge2D>::reinit(num_nodes);
		faces.clear();
		isLeaf.assign(num_nodes, false);
		isBranch.assign(num_nodes, false);
		isInner.assign(num_nodes, false);
		facesfound = false;
		// circumference = 0.0;
	}
	virtual void reinit_flags(int num_nodes)
	{
		isLeaf.assign(num_nodes, false);
		isBranch.assign(num_nodes, false);
		isInner.assign(num_nodes, false);
	}

	void connect(unsigned from, unsigned to, double dx, double dy)
	{
		assert(from < nodes.size());
		assert(to < nodes.size());
		Node2D &fromNode = nodes[from];
		Node2D &toNode = nodes[to];

		Edge2D *edge = getNewEdge();
		edge->setTarget(to);
		edge->setTranslation(dx, dy);
		fromNode.addEdge(edge);

		edge = getNewEdge();
		edge->setTarget(from);
		edge->setTranslation(-dx, -dy);
		toNode.addEdge(edge);
		facesfound = false;
	}

	void connect(unsigned from, unsigned to, v2d dist)
	{
		assert(from < nodes.size());
		assert(to < nodes.size());
		Node2D &fromNode = nodes[from];
		Node2D &toNode = nodes[to];

		Edge2D *edge = getNewEdge();
		edge->setTarget(to);
		edge->setTranslation(dist);
		fromNode.addEdge(edge);

		edge = getNewEdge();
		edge->setTarget(from);
		edge->setTranslation(-dist);
		toNode.addEdge(edge);
		facesfound = false;
	}

	/**
	 * @brief Identify leave nodes (I.e. nodes with only one neighbor) and branch nodes (i.e. nodes connecting to one leaf and one other node) as well as inner nodes (i.e. all other nodes that are adjacent to a face of the graph).
	 *
	 */
	void markNodeTypes()
	{
		reinit_flags(nodes.size());
		// std::cerr << "Marking nodes " << std::endl;
		std::vector<unsigned> leafs;
		std::priority_queue<unsigned> tovisit;
		std::vector<unsigned> num_neighbors(nodes.size());
		std::vector<unsigned> num_neighbors_after_leaves(nodes.size());

		for (unsigned i = 0; i < nodes.size(); i++)
		{
			num_neighbors[i] = nodes[i].getAdjacents().size();
			num_neighbors_after_leaves[i] = num_neighbors[i];
			nodes[i].sortEdges();
			/*std::cerr << "Node #" << i << " sorted nb:" << std::endl;
			for (Edge2D* e : nodes[i].getAdjacents()) {
				std::cerr << "\t" << e->getTarget();
			}
			std::cerr << std::endl;*/
		}
		// std::cerr << "Found neighbor sizes " << std::endl;

		for (unsigned i = 0; i < nodes.size(); i++)
		{
			if (num_neighbors[i] <= 1)
			{
				isLeaf[i] = true;
				leafs.push_back(i);
				for (Edge2D *e : nodes[i].getAdjacents())
				{
					int targ = e->getTarget();
					tovisit.push(targ);
					num_neighbors_after_leaves[targ]--;
					// std::cerr << targ << ": nb=" << num_neighbors[targ] << std::endl;
				}
			}
		}
		// std::cerr << "Found leaves " << std::endl;

		while (!tovisit.empty())
		{
			unsigned curr = tovisit.top();
			tovisit.pop();

			// No need to visit branches/leafs
			if (isLeaf[curr] || isBranch[curr])
			{
				continue;
			}

			// at least two unvisited neighbors means inner/not yet discovered branch
			if (num_neighbors_after_leaves[curr] > 1)
			{
				continue;
			}

			// We have a branch
			isBranch[curr] = true;
			for (Edge2D *e : nodes[curr].getAdjacents())
			{
				int targ = e->getTarget();
				if (!isLeaf[targ] && !isBranch[targ])
				{
					tovisit.push(targ);
					num_neighbors_after_leaves[targ]--;
				}
			}
		}
		// std::cerr << "Found branches " << std::endl;

		// Mark inner nodes
		for (unsigned i = 0; i < nodes.size(); i++)
		{
			isInner[i] = !isLeaf[i] && !isBranch[i];
		}
		// std::cerr << "Marked nodes" << std::endl;
		// this->dump(std::cout);
		// std::cout << std::endl;
	}

	/**
	 * @brief Calculate the area of a face using the cross product sum.
	 *
	 * Does not consider the sign of the orientation but will correctly divide value by two to account for triangulation having half the area of the added cross product.
	 *
	 * @param face
	 * @return double
	 */
	double areaOfFace(Face2D &face)
	{
		double area = 0.0;
		if (face.corners.size() < 3)
		{
			return area;
		}

		v2d base = face.corners[0];
		v2d curr = base + face.edges[0]->getTranslation();

		// Skip first edge already registered and last edge pointing to base
		for (unsigned e = 1; e < face.edges.size() - 1; e++)
		{
			v2d next = curr + face.edges[e]->getTranslation();

			v2d t1 = curr - base;
			v2d t2 = next - base;

			// determinant is double signed area
			area += (t1.real() * t2.imag() - t2.real() * t1.imag());
			curr = next;
		}

		return std::abs(area) / 2.0;
	}
	
	/**
	 * @brief Calculate the perimeter of a face as the sum of edge lengths
	 * 
	 * @param face The face for which to calculate the perimeter
	 * @return double The full perimeter
	*/
	double perimeterOfFace(Face2D &face)
	{
		double perimeter = 0.0;
		if (face.corners.size() < 2)
		{
			return perimeter;
		}

		// Skip first edge already registered and last edge pointing to base
		for (unsigned e = 0; e < face.edges.size(); e++)
		{
			perimeter += face.edges[e]->getLength();
		}

		return perimeter;
	}

	/**
	 * @brief Identify faces in clockwise orientation.
	 *
	 * This is done starting from an "inner" node of the graph and following the outgoing edges back to the original vertex
	 * while only ever turing "right" (i.e. to the next clockwise orientation) edge of the face.
	 *
	 * This is achieved through sorting edges by their angle in each node and taking the "left" edge relative to the ingoing edge when moving to a node.
	 *
	 * @param force_rebuild
	 */
	void findFaces(bool force_rebuild = false)
	{
		if (force_rebuild || !facesfound)
		{
			markNodeTypes();
			// std::cerr << "Finding faces" << std::endl;
			/**
			 * @brief Array to keep track of which edges are considered as already used
			 *
			 * Each edge in either direction must only pe part of one face.
			 */
			std::vector<bool> edge_used(next_edge, false);

			faces.clear();

			// Find all faces involving the marked inner vertices
			for (unsigned origin = 0; origin < nodes.size(); origin++)
			{
				// std::cerr << "start@ " << origin;
				// These edges were sorted in the markNodeTypes() function
				for (Edge2D *e : nodes[origin].getAdjacents())
				{
					/*
						Only consider vertices with at least one neighbor
					*/
					if (!edge_used[edgeindexmap[e]])
					{
						unsigned curr = origin;
						int edgecount = 0;
						Edge2D *current_edge = e;
						Face2D face;

						// We do not want to consider polygons containing the same vertex in two adjacent points.
						// Also we want to iterate until we reach the original point again
						while (current_edge->getTarget() != origin)
						{
							// std::cerr << "\t" << tar;
							// Append current edge to face and mark as used
							edge_used[edgeindexmap[current_edge]] = true;
							face.corners.push_back(curr);
							face.edges.push_back(current_edge);
							edgecount++;

							/*
							 * Try and find the next edge of the face
							 */
							unsigned tar = current_edge->getTarget();
							Edge2D *next_edge = nodes[tar].findLeftof(curr, current_edge->getTranslation());
							// This skips branches and leaves. We therefore removed this.
							/*while (edge_used[edgeindexmap[next_edge]])
							{
								next_edge = nodes[tar].findLeftof(next_edge->getTarget(), -next_edge->getTranslation());
							}*/

							// update iteration and mark edge
							curr = tar;
							current_edge = next_edge;
						}
						// Add the last edge to close the face polygon
						edge_used[edgeindexmap[current_edge]] = true;
						face.corners.push_back(curr);
						face.edges.push_back(current_edge);
						edgecount++;

						// Determine whether the face is internal or external
						//
						face.isInner = isFaceInner(face);

						// Add found face to list of faces.
						faces.emplace_back(face);
					}
				}
			}
			// std::cerr << "Found faces" << std::endl;
			facesfound = true;
		}
	}

	/**
	 * @brief Find the orientation of a face to determine whether it is the outer circumference or an inner face
	 *
	 * We do this by finding the orientation of the complex hull.
	 * The vertex with the lowest x-coordinate is guaranteed to be on the convex hull, so we try and identify that.
	 *
	 * @param face The face to identify as inner or outer
	 * @return true The polygon is clockwise and therefore an inner face.
	 * @return false The polygon is counter-clockwise so considered outside
	 */
	bool isFaceInner(const Face2D &face)
	{
		// Index of the current best candidate for the convex hull
		unsigned convex_hull_index = 0;
		// Use any vertex as initial candidate for the convex hull
		v2d convex_hull_position = nodes[face.corners[convex_hull_index]].getPosition();
		// List to keep track of positions after removing PBC
		std::vector<v2d> pbc_free_nodes;
		pbc_free_nodes.push_back(convex_hull_position);
		// Next position to be checked
		v2d next_possible_position = convex_hull_position;


		// Check all vertices of the polygon whether they are on the convex_hull
		for (unsigned e = 0; e < face.edges.size() - 1; e++)
		{
			next_possible_position = next_possible_position + face.edges[e]->getTranslation();
			pbc_free_nodes.push_back(next_possible_position);
			if (next_possible_position.real() < convex_hull_position.real() || (next_possible_position.real() == convex_hull_position.real() && next_possible_position.imag() < convex_hull_position.imag()))
			{
				convex_hull_index = e + 1;
				convex_hull_position = next_possible_position;
			}
		}
		

		v2d next_position = pbc_free_nodes[(convex_hull_index + 1) % face.corners.size()];
		v2d previous_position = pbc_free_nodes[(convex_hull_index + face.corners.size() - 1) % face.corners.size()];
		double angular_orientation_on_hull = std::arg((next_position - convex_hull_position) / (convex_hull_position - previous_position));

		// Check correct orientation
		bool is_inner = (angular_orientation_on_hull > 0);
		
		// Check total translation of all edges in face
		v2d total_trans;
		for (unsigned e = 0; e < face.edges.size(); e++)
		{
			total_trans += face.edges[e]->getTranslation();
		}
		// Check that edges close up without pbc, which should be true for inner face
		is_inner = is_inner && (abs(total_trans) <  1e-3);

		// If the orientation is counter-clockwise, i.e., angular_orientation_on_hull>0, then this is an inner.
		return is_inner;
	}

	/**
	 * @brief Find the length of the outline of the film
	 *
	 * This is done by adding the identified film outer perimeter to the perimeter of inner holes.
	 * Inner holes are identified based on their number of vertices.
	 *
	 * @param minimum_hole_corners Number of vertices at and above which a face is considered an internal hole of the film.
	 * @return double
	 */
	double findOutline(unsigned minimum_hole_corners = 3)
	{
		findFaces();
		double circumference = 0.0;

		for (unsigned f = 0; f < faces.size(); f++)
		{
			Face2D &face = faces[f];

			// This face is either the outside perimeter polygon or an inside hole with the required size
			if (!face.isInner || face.corners.size() >= minimum_hole_corners)
			{
				double outer_length = 0.0;
				// Only calculate perimeter if required
				for (Edge2D *e : face.edges)
				{
					outer_length += e->getLength();
				}
				circumference += outer_length;
				// std::cerr << "outer face: (" << (!face.isInner) << ": " << outer_length << std::endl;
			}
			else
			{
				// std::cerr << "Discarded outline of length " << outer_length << std::endl;
			}
		}
		// std::cerr << "Found outline" << std::endl;
		return circumference;
	}

	/**
	 * @brief Function to count how many holes of a certain size can be found in the film
	 *
	 * @param minimum_hole_corners Number of vertices at and above which a face is considered an internal hole of the film.
	 * @return std::pair<std::vector<size_t>, size_t> The list of how often a hole with that number of vertices (the index) occurs and the total number of holes as second entry.
	 */
	std::pair<std::vector<size_t>, size_t> findFaceStatistics(unsigned minimum_hole_corners = 3)
	{
		findFaces();

		size_t max_vertex_count = 2 * nodes.size();
		std::vector<size_t> polygon_count(max_vertex_count + 1, 0);
		size_t num_total_holes = 0;

		for (unsigned f = 0; f < faces.size(); f++)
		{
			Face2D &face = faces[f];

			// This face is an inside face. Removed the check for holes at this point as it can be applied in post-processing
			if (face.isInner) // && face.corners.size() >= minimum_hole_corners)
			{
				if (face.corners.size() <= max_vertex_count)
				{
					polygon_count[face.corners.size()]++;
				}
				// Only count holes in number of total holes
				if (face.corners.size() >= minimum_hole_corners)
				{
					num_total_holes++;
				}
			}
		}
		// std::cerr << "Found outline" << std::endl;
		return {polygon_count, num_total_holes};
	}

	/**
	 * @brief Function to find statistics of holes within the graph
	 *
	 * @param minimum_hole_corners Number of vertices at and above which a face is considered an internal hole of the film.
	 * @return std::vector<std::pair<size_t,std::pair<double,double>> List of (vertices,(perimeter,area)) tuples for all holes
	 */
	std::vector<std::pair<size_t,std::pair<double,double>>> findHoleStatistics(unsigned minimum_hole_corners = 3)
	{
		findFaces();

		size_t max_vertex_count = 2 * nodes.size();
		std::vector<std::pair<size_t,std::pair<double,double>>> result;
		size_t num_total_holes = 0;

		for (unsigned f = 0; f < faces.size(); f++)
		{
			Face2D &face = faces[f];

			// This face is an inside face. Removed the check for holes at this point as it can be applied in post-processing
			if (face.isInner) // && face.corners.size() >= minimum_hole_corners)
			{
				// Only count holes in number of total holes
				if (face.corners.size() >= minimum_hole_corners)
				{
					result.push_back({face.corners.size(),{perimeterOfFace(face), areaOfFace(face)}});
				}
			}
		}
		return result;
	}

	/**
	 * @brief Function to determine the area of the graph excluding the area of holes.
	 *
	 * @param minimum_hole_corners Number of vertices at and above which a face is considered an internal hole of the film.
	 * @return double The total film area
	 */
	double findArea(unsigned minimum_hole_corners = 3)
	{
		findFaces();
		double area = 0.0;

		for (unsigned f = 0; f < faces.size(); f++)
		{
			Face2D &face = faces[f];
			// Inside polygon below hole size constraint
			if (face.isInner && face.corners.size() < minimum_hole_corners)
			{
				area += areaOfFace(face);
			}
		}
		return area;
	}

	void dumpNodes(std::ostream &output)
	{
		output << nodes.size() << std::endl;
		for(auto point:nodes){
			output << point.getPosition().real() << "\t" << point.getPosition().imag() << std::endl;
		}
	}

	void dump(std::ostream &output, size_t minimum_hole_corners = 6)
	{
		double maxx = std::numeric_limits<double>::min(), maxy = std::numeric_limits<double>::min();
		double minx = std::numeric_limits<double>::max(), miny = std::numeric_limits<double>::max();
		for (Node2D &n : nodes)
		{
			v2d pos = n.getPosition();
			maxx = std::max(maxx, pos.real());
			maxy = std::max(maxy, pos.imag());
			minx = std::min(minx, pos.real());
			miny = std::min(miny, pos.imag());
		}
		double smallwidth = maxx - minx;
		double smallheight = maxy - miny;
		double width = 800;
		double height = 600;
		double padding = 30;

		output << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << width + 2 * padding << "\" height=\"" << height + 2 * padding << "\">" << std::endl;
		/*output << "<!--g = group; translate(dx dy) : verschieben; "
		<< "rotate(a x y) : a: winkel grad cw, x, y : drehmittelpunkt; "
		<< "scale: zoom, minus = spiegeln(hier : y - koord.fixen)-->";*/
		// output << "<g transform = \"translate(" << -minx << " " << -miny << ") rotate(12.3 0 0) scale(1,-1)\">" << std::endl;

		double wscale = (width / smallwidth);
		double hscale = (height / smallheight);

		double lower_left = padding;
		double lower_bot = padding;

		double scale = std::min(smallwidth, smallheight);
		double pointr = scale / 90.0;

		output << "<g transform = \"translate(" << (lower_left - minx * wscale) << " "
			   << (lower_bot - miny * hscale) << ") scale(" << wscale
			   << ", " << hscale << ")\">" << std::endl;

		/*
		<line y1 = "-50"  x1 = "0" y2 = "400" x2 = "0" stroke = "#808080" / >
		<circle cx = "23" cy = "13" r = "42" fill = "none" stroke = "red" / >
		<polygon fill = "none" stroke = "green" stroke - width = "5"
		points = "100,100 200,100 180,140" / >
		<polyline fill = "lightgreen" stroke = "green" stroke - width = "5"
		points = "100,50 200,50 180,90 220,80" / >
		<ellipse cx = "50" cy = "150" rx = "80" ry = "50" fill = "yellow"
		stroke = "black" stroke - width = "3" stroke - dasharray = "5,10,15,20" / >
		<rect x = "50" y = "150" width = "80" height = "50" fill = "none" stroke = "red" / >
		<rect x = "60" y = "160" width = "80" height = "50" rx = "10" ry = "20"
		fill = "none" stroke = "blue" / >
		<path d = "M-50 -50 L50 -50 L50 50 L-50 50 Z" fill = "none" stroke = "black" / >*/

		findFaces();

		for (Face2D &f : faces)
		{
			if (f.isInner && f.corners.size() < minimum_hole_corners)
			{
				output << "<polyline fill =\"lightgreen\" fill-opacity=\"0.5\" stroke=\"green\" stroke-width=\"" << pointr / 2.0 << "\" "
					   << "points =\"";
				v2d base = nodes[f.corners[0]].getPosition();
				for (Edge2D *e : f.edges)
				{
					base += e->getTranslation();
					output << base.real() << "," << base.imag() << " ";
				}
				output << "\" />" << std::endl;
			}
			else
			{
				output << "<polyline fill =\"lightblue\" fill-opacity=\"0.5\" stroke=\"blue\" stroke-width=\"" << pointr / 2.0 << "\" "
					   << "points =\"";
				v2d base = nodes[f.corners[0]].getPosition();
				for (Edge2D *e : f.edges)
				{
					base += e->getTranslation();
					output << base.real() << "," << base.imag() << " ";
				}
				output << "\" />" << std::endl;
			}
		}

		for (Node2D &n : nodes)
		{
			v2d pos = n.getPosition();
			for (Edge2D *e : n.getAdjacents())
			{
				// Node2D &other = nodes[e->getTarget()];
				v2d opos = pos + e->getTranslation();
				output << "<line x1=\"" << pos.real() << "\" y1 =\"" << pos.imag() << "\" x2 =\"" << opos.real() << "\" y2 =\"" << opos.imag() << "\" stroke=\"black\" stroke-width=\"" << pointr / 2.0 << "\"/>" << std::endl;
			}
		}

		for (unsigned i = 0; i < nodes.size(); i++)
		{
			Node2D &n = nodes[i];
			v2d pos = n.getPosition();
			output << "<circle cx=\"" << pos.real() << "\" cy=\"" << pos.imag() << "\" r=\"" << pointr << "\" fill=\"";
			if (isLeaf[i])
			{
				output << "green";
			}
			else if (isBranch[i])
			{
				output << "blue";
			}
			else
			{
				output << "black";
			}
			output << "\" />" << std::endl;
		}
		output << "</g>" << std::endl;
		output << "</svg>";
	}
};

#endif
