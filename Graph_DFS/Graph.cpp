#include "Graph.h"

Graph::Graph(int V) : V(V), adj(V) 
{
    // Constructor
    // Initialize the adjacency list with the number of vertices
    // Initialize the adjacency list with an empty list for each vertex
} 

/**
 * Adds an undirected edge between two vertices in the graph.
 *
 * This function adds an edge between vertex `v` and vertex `w` in the graph.
 * It updates the adjacency list of both vertices to include each other.
 *
 * @param v The first vertex of the edge.
 * @param w The second vertex of the edge.
 */
void Graph::addEdge(int v, int w) 
{
    adj[v].push_back(w); // Add w to the adjacency list of v.
    adj[w].push_back(v); // Add v to the adjacency list of w since the graph is undirected.
}

int Graph::getNumberOfVertices() const 
{
    return V; // Return the total number of vertices in the graph
}
 
const std::vector<std::list<int>>& Graph::getAdjacencyList() const 
{
    return adj; // Return the adjacency list
}