#ifndef GRAPHREADER_H
#define GRAPHREADER_H

#include "Graph.h"
#include <string>

class GraphReader {
public:
    static Graph readGraphFromFile(const std::string& filename);
};

#endif
