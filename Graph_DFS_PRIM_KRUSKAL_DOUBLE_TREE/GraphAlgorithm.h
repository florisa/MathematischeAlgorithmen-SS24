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
    static std::pair<std::vector<int>, double> doubleTreeTSP(const Graph& graph); 
    static MSTResult primMST(const Graph& graph); 
    static MSTResult kruskalMST(const Graph& graph); 
    static int countConnectedComponents(const Graph& graph);
    static void dfs(int v, std::vector<bool>& visited, const Graph& graph);
};