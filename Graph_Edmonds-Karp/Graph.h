#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

struct Edge {
    int target;
    double weight;
};

class Graph {
public:
    Graph(int numVertices);
    void addEdge(int u, int v, double weight, bool directed = true);
    const std::vector<Edge>& getEdges(int u) const;
    int numVertices() const;

private:
    int numVertices_;
    std::vector<std::vector<Edge>> adjacencyList_;
};

#endif // GRAPH_H
