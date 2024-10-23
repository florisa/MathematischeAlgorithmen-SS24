#ifndef EDGE_H
#define EDGE_H

#include <iostream>

class Edge {
public:
    Edge(int n1, int n2, double w);
    int node1;
    int node2;
    double weight;
};

#endif // EDGE_H