#pragma once

#include "utils.hpp"
#include "graph.hpp"

namespace sspp {
namespace decomp {
TreeDecomposition Treedecomp(const Graph& graph, double time);
} // namespace decomp
} // namespace sspp