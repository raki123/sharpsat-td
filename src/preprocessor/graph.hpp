#pragma once

#include <unordered_set>
#include "utils.hpp"

namespace sspp {

typedef pair<int, int> Edge;

class Graph {
public:
  Graph() : n_(0), m_(0) {};
  Graph(int n) : n_(n), m_(0), adj_list_(n) {};
  Graph(std::vector<Edge> edges);
  Graph(int vars, const vector<vector<Lit>>& clauses);
  Graph(const Graph& other) : n_(other.n_), m_(other.m_), adj_list_(other.adj_list_) {};
  Graph& operator=(const Graph& other);
  ~Graph() {};


  int AddVertex();
  void AddEdge(int v, int u);
  void AddEdge(Edge e);
  void AddEdges(const std::vector<Edge>& edges);

  void RemoveEdge(int v, int u);
  void RemoveEdgesBetween(int v, const std::vector<int>& vs);
  
  int n() const;
  int m() const;
  bool HasEdge(int v, int u) const;
  bool HasEdge(Edge e) const;
  std::vector<Edge> Edges() const;

  
  const std::unordered_set<int>& Neighbors(int v) const;

private:
  int n_, m_;
  std::vector<std::unordered_set<int> > adj_list_;
};

} // namespace sspp