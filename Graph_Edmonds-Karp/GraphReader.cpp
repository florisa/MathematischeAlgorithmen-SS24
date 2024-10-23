#include "GraphReader.h"
#include <fstream>
#include <stdexcept>

Graph GraphReader::readGraphFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    int numVertices;
    file >> numVertices;
    Graph graph(numVertices);

    int u, v;
    double weight;
    while (file >> u >> v >> weight) {
        graph.addEdge(u, v, weight);
    }

    file.close();
    return graph;
}
