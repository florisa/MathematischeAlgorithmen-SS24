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

    try 
    {
        Graph graph = GraphReader::readGraphFromFile(filename); // Read the graph from the specified file
        int componentCount = GraphAlgorithm::countConnectedComponents(graph); // Count the number of connected components in the graph
        std::cout << "Number of connected components: " << componentCount << std::endl; // Output the number of connected components
    } catch (const std::exception& e) 
    {
        std::cerr << "Error: " << e.what() << std::endl; // Output any errors that occurred during the execution
        return 1;
    }

    return 0;
}