#pragma once
#include "Graph.h"
#include "UnionFind.h"  

/**
 * @struct MSTResult
 * @brief Represents the result of a minimum spanning tree algorithm.
 *
 * This class provides functionality to store the result of a minimum spanning tree algorithm.
 */
struct MSTResult 
{
    std::vector<std::pair<int, int>> edges; // Stores pairs of vertices that form the MST
    double totalWeight = 0.0; // Total weight of all the edges in the MST

    // Method to add edges to the MST
    void addMSTEdge(int u, int v, double weight) 
    {
        edges.emplace_back(u, v); // Store the edge as part of the MST result
        totalWeight += weight; // Increment the total weight of the MST
    }

    double getTotalWeight() const 
    {
        return totalWeight; // Getter for the total weight of the MST
    }

    const std::vector<std::pair<int, int>>& getEdges() const 
    {
        return edges; // Getter for the edges list of the MST
    }
};

/**
 * @class GraphAlgorithm
 * @brief Provides implementations of graph algorithms.
 *
 * This class provides implementations of various graph algorithms such as Prim's MST, Kruskal's MST, 
 * counting connected components, and solving the Traveling Salesman Problem (TSP) using the Nearest Neighbor and Double Tree algorithms.
 */

class GraphAlgorithm {
public:
    static std::pair<std::vector<int>, double> nearestNeighborTSP(const Graph& graph);
    static std::pair<std::vector<int>, double> doubleTreeTSP(const Graph& graph); 
    static MSTResult primMST(const Graph& graph); 
    static MSTResult kruskalMST(const Graph& graph); 
    static int countConnectedComponents(const Graph& graph);
    static void dfs(int v, std::vector<bool>& visited, const Graph& graph);

private:
    static Graph doubleEdges(const Graph& graph);
    static std::vector<int> findEulerianTour(const Graph& graph);
};