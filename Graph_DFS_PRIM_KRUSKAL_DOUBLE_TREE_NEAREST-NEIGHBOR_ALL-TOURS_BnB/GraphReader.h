#pragma once

#include "Graph.h"
#include <string>

/**
 * @brief The GraphReader class is responsible for reading a graph from a file.
 */
class GraphReader 
{
public:
    /**
     * @brief Reads a graph from a file.
     * 
     * This function reads a graph from the specified file and returns the graph object.
     * 
     * @param filename The name of the file to read the graph from.
     * @return The graph object read from the file.
     */
    static Graph readGraphFromFile(const std::string &filename);
};