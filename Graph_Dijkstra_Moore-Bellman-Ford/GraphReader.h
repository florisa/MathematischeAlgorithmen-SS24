#ifndef GRAPH_READER_H
#define GRAPH_READER_H

#include "Graph.h"
#include <string>

class GraphReader {
public:
    static Graph readGraphFromFile(const std::string& filename, bool directed);
};

#endif // GRAPH_READER_H
