#include "GraphAlgorithm.h"

/**
 * Counts the number of connected components in the given graph.
 *
 * @param graph The graph to count the connected components in.
 * @return The number of connected components in the graph.
 */
int GraphAlgorithm::countConnectedComponents(const Graph& graph) 
{
    int count = 0;
    std::vector<bool> visited(graph.getNumberOfVertices(), false); // Creates a boolean vector with the size of the number of vertices in the graph to keep track of visited vertices

    for (int i = 0; i < graph.getNumberOfVertices(); i++) 
    {
        if (!visited[i]) 
        {
            dfs(i, visited, graph); // Perform a depth-first search (DFS) traversal of the graph starting from vertex i
            count++; 
        }
    }
    return count;
}

void GraphAlgorithm::dfs(int v, std::vector<bool>& visited, const Graph& graph) 
{
    visited[v] = true; // Mark the current vertex as visited

    for (int neighbor : graph.getAdjacencyList()[v]) // Iterate over the neighbors of vertex v (retrieved from the adjacency list of the graph)
    {
        if (!visited[neighbor]) 
        {
            dfs(neighbor, visited, graph); // Recursively perform a DFS traversal starting from the neighbor
        }
    }
}