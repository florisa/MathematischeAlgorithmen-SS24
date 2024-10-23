#include "Graph.h"
#include "GraphReader.h"
#include "GraphAlgorithm.h"
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

void printDistances(const std::vector<double>& dist) {
    for (size_t i = 0; i < dist.size(); ++i) {
        std::cout << "Vertex " << i << ": " << dist[i] << std::endl;
    }
}

void printPath(int source, int target, const std::vector<int>& pred) {
    std::vector<int> path;
    for (int at = target; at != -1; at = pred[at]) {
        path.push_back(at);
    }
    std::reverse(path.begin(), path.end());
    std::cout << "Path from " << source << " to " << target << ": ";
    for (int vertex : path) {
        std::cout << vertex << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <graph file path> <directed (0 or 1)>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    bool directed = std::stoi(argv[2]) != 0;

    try {
        Graph graph = GraphReader::readGraphFromFile(filename, directed);
        std::cout << "Select an algorithm:\n1. Dijkstra\n2. Moore-Bellman-Ford\n";
        int choice;
        std::cin >> choice;

        std::cout << "Enter the source vertex: ";
        int source;
        std::cin >> source;

        if (choice == 1) {
            auto [distances, predecessors] = GraphAlgorithm::dijkstra(graph, source);
            std::cout << "Distances from source " << source << ":\n";
            printDistances(distances);

            std::cout << "Enter the target vertex to see the path: ";
            int target;
            std::cin >> target;
            printPath(source, target, predecessors);
        } else if (choice == 2) {
            try {
                auto [distances, predecessors] = GraphAlgorithm::mooreBellmanFord(graph, source);
                std::cout << "Distances from source " << source << ":\n";
                printDistances(distances);

                std::cout << "Enter the target vertex to see the path: ";
                int target;
                std::cin >> target;
                printPath(source, target, predecessors);
            } catch (const std::runtime_error& e) {
                std::cerr << "Error: " << e.what() << std::endl;
            }
        } else {
            std::cerr << "Invalid option selected." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
