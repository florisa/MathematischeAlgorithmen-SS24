#pragma once
#include <vector>

class UnionFind 
{
private:
    std::vector<int> parent; // Parent of each element
    std::vector<int> rank; // Rank of each element

public:
    UnionFind(int size); // Constructor to initialize the UnionFind data structure
    int find(int p); 
    void unite(int p, int q); 
};