#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <utility>

struct Edge 
{
    int target; // The target vertex
    double weight;
};

class Graph 
{
private:
    int numVertices_; 
    bool directed_;     
    std::vector<std::vector<Edge>> adjacencyList_; 

public:
    Graph(int numVertices, bool directed);
    void addEdge(int u, int v, double weight);
    const std::vector<Edge>& getEdges(int u) const;
    int numVertices() const;
};

#endif // GRAPH_H
