#include "Graph.h"

Graph::Graph(int V) : V(V), adj(V) 
{
    // Constructor initializes the adjacency list with an empty list for each vertex
}

void Graph::addEdge(int v, int w, double weight) 
{
    adj[v].push_back({w, weight}); // Add w and weight to the adjacency list of v
    adj[w].push_back({v, weight}); // Add v and weight to the adjacency list of w since the graph is undirected
}

int Graph::getNumberOfVertices() const 
{
    return V; // Return the total number of vertices in the graph
}

const std::vector<std::list<std::pair<int, double>>>& Graph::getAdjacencyList() const 
{
    return adj; // Return the adjacency list
}
