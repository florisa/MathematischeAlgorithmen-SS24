#include "GraphReader.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>

/**
 * Reads a graph from a file.
 * 
 * @param filename The name of the file to read the graph from.
 * @param directed A boolean indicating whether the graph is directed or not.
 * @return The graph read from the file.
 * @throws std::runtime_error if the file cannot be opened.
 */
Graph GraphReader::readGraphFromFile(const std::string& filename, bool directed) 
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        throw std::runtime_error("Could not open file");
    }

    int numVertices;
    file >> numVertices;

    Graph graph(numVertices, directed);
    int u, v;
    double weight;
    while (file >> u >> v >> weight) 
    {
        graph.addEdge(u, v, weight);
    }

    file.close();

    // Debug *****************************************************************************************************
    /* std::cout << "Graph read from file:\n";
    for (int i = 0; i < numVertices; ++i) 
    {
        for (const auto& edge : graph.getEdges(i)) 
        {
            std::cout << "Edge from " << i << " to " << edge.target << " with weight " << edge.weight << std::endl;
        }
    } */

    return graph;
}
