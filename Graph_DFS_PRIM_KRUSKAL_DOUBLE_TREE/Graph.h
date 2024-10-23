#pragma once
#include <vector> 
#include <list> 
#include <utility> // For std::pair

/**
 * @class Graph
 * @brief Represents a weighted graph data structure.
 *
 * This class provides functionality to create and manipulate a weighted graph.
 */
class Graph 
{
private:
    int V; // Number of vertices
    std::vector<std::list<std::pair<int, double>>> adj; // Adjacency list with weights

public:
    /** Constructor to initialize the graph with a specified number of vertices */
    Graph(int V);

    /** Adds an undirected edge between two vertices with a weight */
    void addEdge(int v, int w, double weight);


    /** Returns the number of vertices */
    int getNumberOfVertices() const;

    /** Provides access to the adjacency list */
    const std::vector<std::list<std::pair<int, double>>>& getAdjacencyList() const;

    /** Helper Methods for Double Tree TSP*/
    Graph doubleEdges() const;
    std::vector<int> findEulerianTour() const;
};