#include "GraphAlgorithm.h"
#include <algorithm> 
#include <queue>     
#include <set>  


// Prim's algorithm implementation to compute the MST of a graph
MSTResult GraphAlgorithm::primMST(const Graph& graph) 
{
    int V = graph.getNumberOfVertices(); // Get the number of vertices in the graph
    std::vector<bool> inMST(V, false); // Vector to keep track of vertices in the MST, initially all false
   
    /* Declare a priority queue 'pq' using a container the std::priority_queue, 
       with the type std::pair<double, int> (weight, vertex) and a custom comparator 
       std::greater<> to get the smallest weight edge on top (smallest elements have the highest priority).
       //https://stackoverflow.com/questions/42100055/stdpriority-queue-with-stdpairint-int
    */

    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;
    MSTResult mstResult;

    // Start with vertex 0 
    inMST[0] = true;
    const auto& adjList = graph.getAdjacencyList(); //retrieves the adjacency list of the graph
    for (const auto& edge : adjList[0]) 
    {  // loop through all the edges of vertex 0
        pq.emplace(edge.second, edge.first); // For each edge, it adds a pair (the edge's weight and the vertex it is connected to) to the priority queue pq
    }

    while (!pq.empty()) 
    { 
        auto [weight, u] = pq.top(); // Get the edge with the smallest weight
        pq.pop(); // Remove the edge from the priority queue

        if (!inMST[u]) // If the vertex is not in the MST
        { 
            inMST[u] = true; // Add the vertex to the MST
            mstResult.addMSTEdge(pq.top().second, u, weight); // Add the edge to the MST result

            for (const auto& edge : adjList[u]) // Iterate over all edges of the vertex u
            { 
                if (!inMST[edge.first]) // If the vertex is not in the MST
                { 
                    pq.emplace(edge.second, edge.first); // Add the edge to the priority queue
                }
            }
        }
    }

    return mstResult;
}

// Kruskal's algorithm using a union-find data structure
MSTResult GraphAlgorithm::kruskalMST(const Graph& graph) 
{
    std::vector<std::tuple<double, int, int>> edges; // weight, start vertex, end vertex
    const auto& adjList = graph.getAdjacencyList(); // Get the adjacency list of the graph

    for (int v = 0; v < adjList.size(); ++v) // Iterate over all vertices
    { 
        for (const auto& edge : adjList[v]) // Loop through all edges of the vertex
        { 
            if (v < edge.first) // Only consider edges once (undirected graph), ensure each edge is only added once!!
            { 
                edges.emplace_back(edge.second, v, edge.first); 
            }
        }
    }

    // Sort edges by ascending weight
    std::sort(edges.begin(), edges.end());

    UnionFind uf(graph.getNumberOfVertices()); // Initialize a union-find data structure with the number of vertices
    MSTResult mstResult;

    for (const auto& [weight, u, v] : edges) // Iterate over sorted edges
    { 
        if (uf.find(u) != uf.find(v)) // Check if adding the edge creates a cycle
        { 
            uf.unite(u, v); // If not, add the edge to the MST
            mstResult.addMSTEdge(u, v, weight); // Add the edge to the MST result
            if (mstResult.getEdges().size() == graph.getNumberOfVertices() - 1) break; // Stop when MST is complete
        }
    }

    return mstResult;
}

// Count the number of connected components in the graph using DFS
int GraphAlgorithm::countConnectedComponents(const Graph& graph) 
{
    int V = graph.getNumberOfVertices();
    std::vector<bool> visited(V, false);
    int count = 0;

    for (int i = 0; i < V; i++) 
    {
        if (!visited[i]) 
        {
            dfs(i, visited, graph);
            count++;
        }
    }

    return count;
}

// Depth-first search to explore the graph
void GraphAlgorithm::dfs(int v, std::vector<bool>& visited, const Graph& graph) 
{
    visited[v] = true;
    for (const auto& edge : graph.getAdjacencyList()[v]) 
    {
        if (!visited[edge.first]) 
        {
            dfs(edge.first, visited, graph);
        }
    }
}