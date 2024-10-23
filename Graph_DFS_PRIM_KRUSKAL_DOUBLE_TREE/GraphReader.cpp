#include "GraphReader.h"
#include <fstream>
#include <sstream>

Graph GraphReader::readGraphFromFile(const std::string& filename) {
    std::ifstream file(filename);
    int V;
    std::string line;

    if (!file) {
        throw std::runtime_error("File could not be opened.");
    }

    std::getline(file, line);
    std::istringstream iss(line);
    iss >> V;  // Extract the number of vertices from the line

    Graph graph(V);  // Create a graph with the specified number of vertices

    while (std::getline(file, line)) {
        std::istringstream edgeStream(line);
        int v, w;
        double weight;
        if (edgeStream >> v >> w >> weight) {
            graph.addEdge(v, w, weight);  // Read weighted edge
        } else {
            edgeStream.clear();  // Clear the failed state after attempting to read weight
            edgeStream.seekg(0);  // Reset stream position to re-read the line
            edgeStream >> v >> w;
            graph.addEdge(v, w, 1.0);  // Assume default weight of 1 for unweighted edges
        }
    }

    return graph;
}
