#ifndef GRAPH_ALGORITHM_H
#define GRAPH_ALGORITHM_H

#include "Graph.h"
#include <vector>
#include <utility>

class GraphAlgorithm {
public:
    static std::pair<std::vector<double>, std::vector<int>> dijkstra(const Graph& graph, int source);
    static std::pair<std::vector<double>, std::vector<int>> mooreBellmanFord(const Graph& graph, int source);
};

#endif // GRAPH_ALGORITHM_H
