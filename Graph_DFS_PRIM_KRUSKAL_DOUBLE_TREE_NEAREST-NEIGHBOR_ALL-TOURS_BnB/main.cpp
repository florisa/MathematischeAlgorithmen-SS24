#include "Graph.h"
#include "GraphReader.h"
#include "GraphAlgorithm.h"
#include <iostream>
#include <string>
#include <chrono>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <graph file path>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    try {
        Graph graph = GraphReader::readGraphFromFile(filename);
        std::cout << "Select an operation:\n1. Solve TSP using Branch and Bound\n2. Solve TSP using All Tours\n";
        int choice;
        std::cin >> choice;

        switch (choice) {
            case 1: {
                auto start = std::chrono::high_resolution_clock::now(); // Start timing
                double totalWeight = GraphAlgorithm::solveTSPBranchAndBound(graph);
                auto end = std::chrono::high_resolution_clock::now(); // End timing
                std::chrono::duration<double> duration = end - start;
                std::cout << "Total weight of the TSP tour (Branch-and-Bound): " << totalWeight << std::endl;
                std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
                break;
            }
            case 2: {
                auto start = std::chrono::high_resolution_clock::now(); // Start timing
                double totalWeight = GraphAlgorithm::allTours(graph);
                auto end = std::chrono::high_resolution_clock::now(); // End timing
                std::chrono::duration<double> duration = end - start;
                std::cout << "Total weight of the TSP tour (All Tours): " << totalWeight << std::endl;
                std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
                break;
            }
            default:
                std::cout << "Invalid option selected." << std::endl;
                break;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
