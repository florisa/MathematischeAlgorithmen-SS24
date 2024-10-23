#include "Graph.h"
#include "GraphReader.h"
#include "GraphAlgorithm.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <graph file path>" << std::endl;
        return 1;
    }

    try {
        Graph graph = GraphReader::readGraphFromFile(argv[1]);

        std::cout << "Enter the source vertex: ";
        int source;
        std::cin >> source;

        std::cout << "Enter the sink vertex: ";
        int sink;
        std::cin >> sink;

        int maxFlow = GraphAlgorithm::edmondsKarp(graph, source, sink);
        std::cout << "The maximum possible flow is " << maxFlow << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
