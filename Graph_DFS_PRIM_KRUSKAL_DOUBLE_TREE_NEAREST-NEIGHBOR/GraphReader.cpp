#include "GraphReader.h"
#include <fstream>
#include <sstream>

/**
 * Reads a graph from a file.
 * @param filename The name of the file to read the graph from.
 * @return The graph read from the file.
 * @throws std::runtime_error if the file could not be opened.
 */
Graph GraphReader::readGraphFromFile(const std::string& filename) 
{
    std::ifstream file(filename);
    int V;
    std::string line;

    if (!file) 
    {
        throw std::runtime_error("File could not be opened.");
    }

    std::getline(file, line); // Read the first line of the file (Here assumed to contain the number of vertices)
    std::istringstream iss(line); // Create a stream from the first line 
    iss >> V;  // Extract the number of vertices from the line

    Graph graph(V);  // Create a graph with the specified number of vertices

    while (std::getline(file, line)) 
    {
        std::istringstream edgeStream(line); // Create a stream from the line 
        int v, w; 
        double weight;

        if (edgeStream >> v >> w >> weight) // Extract the vertices and weight from the line
        {
            graph.addEdge(v, w, weight);  
        } else // If the line does not contain a weight
        {
            edgeStream.clear();  // Clear the failed state after attempting to read weight
            edgeStream.seekg(0);  // Reset stream position to re-read the line
            edgeStream >> v >> w; // Extract the vertices without weight
            graph.addEdge(v, w, 1.0);  // Assume default weight of 1 for unweighted edges
        }
    }

    return graph;
}
