#include "graph.hpp"

#include "utils.hpp"

#include <queue>

namespace sspp {


Graph::Graph(std::vector<Edge> edges) : n_(0), m_(0) {
  for(auto edge : edges) {
    int highest = max(edge.first, edge.second);
    if(highest > n_) {
      n_ = highest + 1;
      adj_list_.resize(n_ + 1);
      AddEdge(edge);
    }
  }
}

Graph::Graph(int vars, const vector<vector<Lit>>& clauses) : Graph((int)vars+1) {
	for (const auto& clause : clauses) {
		for (int i = 0; i < (int)clause.size(); i++) {
			for (int j = i+1; j < (int)clause.size(); j++) {
				Var v1 = VarOf(clause[i]);
				Var v2 = VarOf(clause[j]);
				assert(v1 >= 1 && v1 <= vars && v2 >= 1 && v2 <= vars);
				if (v1 != v2) {
					AddEdge(v1, v2);
				}
			}
		}
	}
}

Graph& Graph::operator=(const Graph& other) {
  n_ = other.n_;
  m_ = other.m_;
  adj_list_ = other.adj_list_;
  return *this;
}

int Graph::n() const {
  return n_;
}

int Graph::m() const {
  return m_;
}

bool Graph::HasEdge(int v, int u) const {
  return std::find(adj_list_[v].begin(), adj_list_[v].end(), u) != adj_list_[v].end();
}

bool Graph::HasEdge(Edge e) const {
  return HasEdge(e.first, e.second);
}

std::vector<Edge> Graph::Edges() const {
  std::vector<Edge> ret;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_list_[i]) {
      if (a > i) ret.push_back({i, a});
    }
  }
  return ret;
}

const std::unordered_set<int>& Graph::Neighbors(int v) const {
  return adj_list_[v];
}

void Graph::AddEdge(int v, int u) {
  assert(v != u);
  if(HasEdge(u,v)) return;
  m_++;
  adj_list_[v].insert(u);
  adj_list_[u].insert(v);
}

void Graph::AddEdge(Edge e) {
  AddEdge(e.first, e.second);
}

void Graph::AddEdges(const std::vector<Edge>& edges) {
  for (auto& edge : edges) AddEdge(edge);
}

void Graph::RemoveEdge(int v, int u) {
  assert(HasEdge(v, u) && HasEdge(u, v));
  m_--;
  adj_list_[u].erase(v);
  adj_list_[v].erase(u);
}

} // namespace sspp