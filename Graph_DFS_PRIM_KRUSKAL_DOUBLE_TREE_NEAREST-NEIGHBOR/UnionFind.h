#pragma once
#include <vector>

class UnionFind 
{
private:
    std::vector<int> parent; //This vector stores the parent of each element, if parent[i] = i, then i is a root node
    std::vector<int> rank; // This vector is used for union by rank, it keeps track of the depth of the tree

public:
    UnionFind(int size); // Constructor to initialize the UnionFind data structure
    int find(int p); 
    void unite(int p, int q); 
};