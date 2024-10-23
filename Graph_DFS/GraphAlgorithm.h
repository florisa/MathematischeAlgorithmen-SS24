#pragma once
#include "Graph.h"
#include <vector>

/**
 * @class GraphAlgorithm
 * @brief A class that provides algorithms for working with graphs.
 */
class GraphAlgorithm 
{
public:
    /**
     * @brief Counts the number of connected components in a graph.
     * @param graph The graph to count the connected components in.
     * @return The number of connected components in the graph.
     */
    static int countConnectedComponents(const Graph& graph); // This function takes a graph as input and returns the number of connected components in the graph.

private:
    /**
     * @brief Performs a depth-first search (DFS) traversal of the graph.
     * @param v The vertex to start the DFS traversal from.
     * @param visited A vector to keep track of visited vertices.
     * @param graph The graph to perform the DFS traversal on.
     */
    static void dfs(int v, std::vector<bool>& visited, const Graph& graph); // This function performs a depth-first search (DFS) traversal of the graph.
};