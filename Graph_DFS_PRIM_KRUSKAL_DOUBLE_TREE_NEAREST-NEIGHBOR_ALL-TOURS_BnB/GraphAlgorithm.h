#ifndef GRAPH_ALGORITHM_H
#define GRAPH_ALGORITHM_H

#include "Graph.h"
#include <vector>
#include <limits>
#include <queue>
#include <algorithm>
#include <numeric>

class GraphAlgorithm {
public:
    static double solveTSPBranchAndBound(const Graph& graph);
    static double allTours(const Graph& graph);

private:
    struct Node 
    {
        std::vector<int> path;
        double cost;
        double lowerBound;
    };

    struct CompareNode 
    {
        bool operator()(const Node& lhs, const Node& rhs) const 
        {
            return lhs.lowerBound > rhs.lowerBound; 
        }
    };

    static double calculateTourCost(const Graph& graph, const std::vector<int>& tour);
    static double calculateLowerBound(const Graph& graph, const std::vector<int>& path);
};

#endif // GRAPH_ALGORITHM_H
