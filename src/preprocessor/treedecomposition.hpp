#pragma once

#include "graph.hpp"

namespace sspp {


class TreeDecomposition {
 public:
    TreeDecomposition() : bs(0), n(0), width(-1) {};
 	TreeDecomposition(const Graph& graph, double time);
    TreeDecomposition(const TreeDecomposition& other) : bs(other.bs), n(other.n), width(other.width), tree(other.tree), bags(other.bags) {};
    TreeDecomposition& operator=(const TreeDecomposition& other);
    ~TreeDecomposition() {};
 	vector<vector<int>>& Bags();
 	const std::unordered_set<int>& Neighbors(int b) const;
 	int nverts() const;
 	int nbags() const;
 	int Width() const;
 	bool Verify(const Graph& graph) const;
 	bool InBag(int b, int v) const;
 	Graph Chordal() const;
 	int Centroid() const;
    vector<int> GetOrd() const;
 	int bs, n, width;
 private:
 	Graph tree;
 	vector<vector<int>> bags;
    void OdDes(int b, int p, int d, vector<int>& ret) const;
 	int CenDfs(int x, int p, int& cen) const;
 	bool dfs(int x, int v, int p, vector<int>& u) const;
};

}