#pragma once
#include "Graph.h"
#include "UnionFind.h"  

struct MSTResult 
{
    std::vector<std::pair<int, int>> edges; // Stores pairs of vertices that form the MST
    double totalWeight = 0.0; // Total weight of all the edges in the MST

    // Method to add edges to the MST
    void addMSTEdge(int u, int v, double weight) 
    {
        /* push_back() vs emplace_back()
        * push_back() creates a temporary object and then copies it to the vector.
        * emplace_back() constructs the object directly in the vector, which is more efficient/faster.
        * In this case, we are constructing a pair of integers, so emplace_back() is more efficient.
        * https://www.geeksforgeeks.org/push_back-vs-emplace_back-in-cpp-stl-vectors/
        */
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

class GraphAlgorithm 
{
public:
    static MSTResult primMST(const Graph& graph); 
    static MSTResult kruskalMST(const Graph& graph); 
    static int countConnectedComponents(const Graph& graph);
    static void dfs(int v, std::vector<bool>& visited, const Graph& graph);
};