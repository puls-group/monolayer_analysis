#ifndef __GRAPH_H__
#define __GRAPH_H__
#include <cmath>
#include <map>
#include <vector>
//#include <iostream>

class Edge;

template<class EDGECLASS>
class Node;

template <class NODECLASS, class EDGECLASS>
class Graph;

class Edge {
protected:
	unsigned _target;
public:
	Edge() :_target(-1) {}

	Edge(unsigned target) :_target(target) {}

	unsigned getTarget() {
		return _target;
	}

	void setTarget(unsigned target) {
		_target = target;
	}
};

template<class EDGECLASS>
class Node {
protected:
	std::vector<EDGECLASS*> _adjacent;
public:
	Node() {}
	virtual ~Node() {}

	void reset() {
		_adjacent.clear();
	}

	virtual void addEdge(EDGECLASS *edge) {
		_adjacent.push_back(edge);
	}

	virtual std::vector<EDGECLASS*> &getAdjacents() {
		return _adjacent;
	}
};

template <class NODECLASS, class EDGECLASS>
class Graph {
protected:
	std::vector<NODECLASS> nodes;
	std::vector<EDGECLASS*> edges;
	std::map<EDGECLASS*, int> edgeindexmap;
	unsigned next_edge;
public:
	Graph(int num_nodes) :
		nodes(0), edges(0), next_edge(0)
	{
		this->reinit(num_nodes);
	}

	virtual ~Graph() {
		for (EDGECLASS * e : edges) {
			delete e;
		}
	}

	virtual void reinit(int num_nodes) {
		nodes.resize(num_nodes);
		for (unsigned n = 0; n < nodes.size(); n++) {
			nodes[n].reset();
		}
		unsigned new_size = num_nodes*num_nodes;
		unsigned old_size = edges.size();
		if (old_size > new_size) {
			for (unsigned e = new_size; e < edges.size(); e++) {
				delete edges[e];
			}
			edges.resize(new_size);
		}
		else {
			edges.resize(new_size);
			for (unsigned e = old_size; e < new_size; e++) {
				edges[e] = new EDGECLASS();
			}
		}
		next_edge = 0;
		edgeindexmap.clear();
	}

	EDGECLASS * getNewEdge() {
		if (next_edge >= edges.size()) {
			edges.resize(edges.size() + 10);
			for (int i = 0; i < 10; i++) {
				edges[next_edge + i] = new EDGECLASS();
			}
		}
		EDGECLASS *result = edges[next_edge];
		edgeindexmap[result] = next_edge;
		next_edge++;
		return result;
	}

	NODECLASS *getNode(int i) {
		if (i < 0 || unsigned(i) >= nodes.size()) {
			return nullptr;
		}
		return &nodes[i];
	}
};
#endif
