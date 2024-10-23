#pragma once
#include <vector>
#include <list>

/**
 * @class Graph
 * @brief Represents a graph data structure.
 *
 * This class provides functionality to create and manipulate a graph.
 */
class Graph {
private:
    int V; // Number of vertices
    std::vector<std::list<int>> adj; // Adjacency list

public:
    /**
     * @brief Constructs a graph with the specified number of vertices.
     * @param V The number of vertices in the graph.
     */
    Graph(int V);

    /**
     * @brief Adds an edge between two vertices in the graph.
     * @param v The source vertex.
     * @param w The destination vertex.
     */
    void addEdge(int v, int w);

    /**
     * @brief Gets the number of vertices in the graph.
     * @return The number of vertices.
     */
    int getNumberOfVertices() const;

    /**
     * @brief Gets the adjacency list of the graph.
     * @return The adjacency list.
     */
    const std::vector<std::list<int>>& getAdjacencyList() const;
};
/* class Graph {
private:
    int V; // Number of vertices
    std::vector<std::list<int>> adj; // Adjacency list

public:
    Graph(int V); // Constructor
    void addEdge(int v, int w); // Function to add an edge to the graph
    int getNumberOfVertices() const; // Getter for the number of vertices
    const std::vector<std::list<int>>& getAdjacencyList() const; // Accessor for adjacency list
}; */