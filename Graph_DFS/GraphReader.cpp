#include "GraphReader.h"
#include <fstream>
#include <sstream>

/**
 * Reads a graph from a file.
 *
 * This function reads a graph from the specified file and returns the graph object.
 *
 * @param filename The name of the file to read the graph from.
 * @return The graph object read from the file.
 */
Graph GraphReader::readGraphFromFile(const std::string& filename) 
{
    std::ifstream file(filename); 
    int V, v, w; 
    std::string line; 

    if (!file) // Check if the file could be opened
    {
        throw std::runtime_error("File could not be opened."); 
    }

    std::getline(file, line); //reads a single line from an input file and stores it in the variable line
    std::istringstream iss(line); // Create a string stream from the line (parse the content of line)
    iss >> V; // Extract the number of vertices from the line

    Graph graph(V); // Create a graph with the specified number of vertices

    while (std::getline(file, line)) 
    {
        std::istringstream edgeStream(line); // Create a string stream from the line
        edgeStream >> v >> w; // Extract the vertices of the edge (v, the source vertex, and w, the target vertex)
        graph.addEdge(v, w); // Add the edge to the graph
    }

    return graph; 
}


/* Graph GraphReader::readGraphFromFile(const std::string& filename) 
{
    std::ifstream file(filename); 
    int v, w, maxVertex = 0; 
    std::string line; // Read the file line by line

    if (!file) // Check if the file could be opened
    {
        throw std::runtime_error("File could not be opened."); 
    }

    while (std::getline(file, line)) 
    {
        std::istringstream edgeStream(line); // Create a string stream from the line
        edgeStream >> v >> w; // Read the vertices of the edge
        maxVertex = std::max(maxVertex, std::max(v, w)); // Update the maximum vertex
    }

    Graph graph(maxVertex + 1); // Create a graph with the maximum vertex

    file.clear(); // Clear the end-of-file flag
    file.seekg(0, std::ios::beg); // Move the file pointer to the beginning of the file

    while (std::getline(file, line)) // Read the file line by line
    {
        std::istringstream edgeStream(line); // Create a string stream from the line
        edgeStream >> v >> w; // Read the vertices of the edge
        graph.addEdge(v, w); // Add the edge to the graph
    }

    return graph;
} */