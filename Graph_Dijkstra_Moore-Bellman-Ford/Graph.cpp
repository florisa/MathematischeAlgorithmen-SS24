#include "Graph.h"

/**
 * @brief Constructs a Graph object.
 * 
 * This constructor initializes a Graph object with the specified number of vertices and
 * whether the graph is directed or not. It also initializes the adjacency list with the
 * specified number of vertices.
 * 
 * @param numVertices The number of vertices in the graph.
 * @param directed A boolean value indicating whether the graph is directed or not.
 */
Graph::Graph(int numVertices, bool directed) 
    : numVertices_(numVertices), directed_(directed), adjacencyList_(numVertices) // Member initializer list
    {
        // Nothing to do here
    }

void Graph::addEdge(int u, int v, double weight) 
{
    adjacencyList_[u].push_back({v, weight});
    if (!directed_) // If the graph is undirected, add the reverse edge
    {
        adjacencyList_[v].push_back({u, weight}); 
    }
}

const std::vector<Edge>& Graph::getEdges(int u) const 
{
    return adjacencyList_[u];
}

int Graph::numVertices() const 
{
    return numVertices_;
}
