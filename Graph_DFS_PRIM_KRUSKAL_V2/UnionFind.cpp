#include "UnionFind.h"

UnionFind::UnionFind(int size) : parent(size), rank(size, 0) 
{
    for (int i = 0; i < size; ++i) {
        parent[i] = i; // Each element is initially its own parent
    }
}

int UnionFind::find(int p) 
{
    if (parent[p] != p) // Check if the element is its own parent (if an element is its own parent, it means that the element is the representative member (or the “root”) of its set.)
    {
        parent[p] = find(parent[p]); // Path compression (each node in the path from p to the root is updated to point to the root directly, reducing the height of the tree and improving the efficiency of future find operations)
    }
    return parent[p]; // Return the parent of p (the root of the set that p belongs to)
}

void UnionFind::unite(int p, int q) // Unite the sets that contain 'p' and 'q' 
{ 
    int rootP = find(p); // Find the root of the set that 'p' belongs to 
    int rootQ = find(q); // Find the root of the set that 'q' belongs to

    if (rootP != rootQ) // If 'p' and 'q' are already in the same set, do nothing
    {
        if (rank[rootP] < rank[rootQ]) // If the rank of 'p' is less than the rank of 'q', make 'q' the parent of 'p'
        { 
            parent[rootP] = rootQ;
        } else if (rank[rootP] > rank[rootQ]) // If the rank of 'p' is greater than the rank of 'q', do the opposite
        { 
            parent[rootQ] = rootP;
        } else // If the ranks are equal, make one a child of the other and increment the rank
        {
            parent[rootQ] = rootP; 
            rank[rootP]++; 
        }
    }
}