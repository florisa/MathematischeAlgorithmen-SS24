#include "Graph.h"

Graph::Graph(int V) : V(V), adj(V) 
{
    // Constructor initializes the adjacency list with an empty list for each vertex
}

void Graph::addEdge(int v, int w, double weight) // Populate the adjacency list with the vertices and weights of the edges
{
    adj[v].push_back({w, weight}); // Add w and weight to the adjacency list of v
    adj[w].push_back({v, weight}); // Add v and weight to the adjacency list of w since the graph is undirected
}

int Graph::getNumberOfVertices() const 
{
    return V; // Return the total number of vertices in the graph
}

const std::vector<std::list<std::pair<int, double>>>& Graph::getAdjacencyList() const // Vector where each elment is a list of pairs (vertex, weight), this stores the adjacency list
{
    return adj; 
}

double Graph::getEdgeWeight(int u, int v) const // Retrieve the weight of the edge between two vertices
{
    for (const auto& [vertex, weight] : adj[u]) //Query the adjacency list of vertex u
    {
        if (vertex == v) // If the vertex is found in the adjacency list
        {
            return weight;
        }
    }
    throw std::out_of_range("Edge does not exist");
}