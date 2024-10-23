#ifndef GRAPH_ALGORITHM_H
#define GRAPH_ALGORITHM_H

#include "Graph.h"

class GraphAlgorithm {
public:
    static int edmondsKarp(const Graph& graph, int source, int sink);
};

#endif // GRAPH_ALGORITHM_H
