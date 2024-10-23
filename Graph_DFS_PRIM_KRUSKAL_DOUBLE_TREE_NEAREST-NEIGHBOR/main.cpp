#include "Graph.h"
#include "GraphReader.h"
#include "GraphAlgorithm.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) 
{
    if (argc != 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <graph file path>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    try {
        Graph graph = GraphReader::readGraphFromFile(filename);
        std::cout << "Select an operation:\n1. Prim's MST\n2. Kruskal's MST\n3. Count Connected Components\n4. Double Tree TSP\n5. Nearest Neighbor TSP\n";
        int choice;
        std::cin >> choice;

        switch (choice) 
        {
            case 1: {
                MSTResult result = GraphAlgorithm::primMST(graph);
                std::cout << "Total weight of MST (Prim's): " << result.getTotalWeight() << std::endl;
                break;
            }
            case 2: {
                MSTResult result = GraphAlgorithm::kruskalMST(graph);
                std::cout << "Total weight of MST (Kruskal's): " << result.getTotalWeight() << std::endl;
                break;
            }
            case 3: {
                int components = GraphAlgorithm::countConnectedComponents(graph);
                std::cout << "Number of connected components: " << components << std::endl;
                break;
            }
            case 4: {
                auto [tspTour, totalWeight] = GraphAlgorithm::doubleTreeTSP(graph);
                std::cout << "Total weight of the TSP tour (Double Tree Algorithm): " << totalWeight << std::endl;
                break;
            }
            case 5: {
                auto [tspTour, totalWeight] = GraphAlgorithm::nearestNeighborTSP(graph);
                std::cout << "Total weight of the TSP tour (Nearest Neighbor Algorithm): " << totalWeight << std::endl;
                break;
            }
            default:
                std::cout << "Invalid option selected." << std::endl;
                break;
        }
    } catch (const std::exception& e) 
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
