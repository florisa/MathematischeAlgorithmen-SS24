#include "Graph.h"

Graph::Graph(int numVertices)
    : numVertices_(numVertices), adjacencyList_(numVertices) {}

void Graph::addEdge(int u, int v, double weight, bool directed) {
    adjacencyList_[u].push_back({v, weight});
    if (!directed) {
        adjacencyList_[v].push_back({u, weight});
    }
}

const std::vector<Edge>& Graph::getEdges(int u) const {
    return adjacencyList_[u];
}

int Graph::numVertices() const {
    return numVertices_;
}
