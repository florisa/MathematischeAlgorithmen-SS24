#include "GraphAlgorithm.h"
#include <algorithm>
#include <queue>
#include <set>
#include <stack>
#include <unordered_set>
#include <vector>
#include <list>
#include <unordered_map>
#include <limits> 

// Helper function to create a doubled graph
Graph doubleEdges(const Graph& graph) {
    Graph doubledGraph(graph.getNumberOfVertices());
    const auto& adj = graph.getAdjacencyList();

    for (int v = 0; v < graph.getNumberOfVertices(); ++v) {
        for (const auto& [w, weight] : adj[v]) {
            doubledGraph.addEdge(v, w, weight);
            doubledGraph.addEdge(w, v, weight);  // Add both directions
        }
    }

    return doubledGraph;
}

// Helper function to find an Eulerian tour
std::vector<int> findEulerianTour(const Graph& graph) {
    std::vector<int> tour;
    std::stack<int> stack;
    std::unordered_map<int, std::multiset<int>> adjList;

    const auto& adj = graph.getAdjacencyList();
    for (int v = 0; v < graph.getNumberOfVertices(); ++v) {
        for (const auto& [w, weight] : adj[v]) {
            adjList[v].insert(w);
        }
    }

    int start = 0;  // Starting vertex
    stack.push(start);

    while (!stack.empty()) {
        int v = stack.top();
        if (adjList[v].empty()) {
            tour.push_back(v);
            stack.pop();
        } else {
            int u = *adjList[v].begin();
            stack.push(u);
            adjList[v].erase(adjList[v].begin());
            adjList[u].erase(adjList[u].find(v));
        }
    }

    return tour;
}

// Double Tree Algorithm for TSP approximation
std::pair<std::vector<int>, double> GraphAlgorithm::doubleTreeTSP(const Graph& graph) {
    // Step 1: Construct the MST using Prim's algorithm
    MSTResult mstResult = kruskalMST(graph);
    Graph mstGraph(graph.getNumberOfVertices());
    for (const auto& edge : mstResult.edges) {
        int u = edge.first;
        int v = edge.second;
        // Find the weight of the edge (u, v) in the original graph
        double weight = -1;
        for (const auto& [w, edgeWeight] : graph.getAdjacencyList()[u]) {
            if (w == v) {
                weight = edgeWeight;
                break;
            }
        }
        mstGraph.addEdge(u, v, weight);
    }

    // Step 2: Double the edges of the MST
    Graph doubledMST = doubleEdges(mstGraph);

    // Step 3: Compute an Eulerian tour of the doubled MST
    std::vector<int> eulerianTour = findEulerianTour(doubledMST);

    // Step 4: Transform the Eulerian tour into a Hamiltonian cycle
    std::vector<int> hamiltonianCycle;
    std::unordered_set<int> visited;
    double totalWeight = 0.0;

    int prevVertex = -1;
    for (int vertex : eulerianTour) {
        if (visited.find(vertex) == visited.end()) {
            if (prevVertex != -1) {
                // Find the weight of the edge (prevVertex, vertex) in the original graph
                for (const auto& [w, weight] : graph.getAdjacencyList()[prevVertex]) {
                    if (w == vertex) {
                        totalWeight += weight;
                        break;
                    }
                }
            }
            visited.insert(vertex);
            hamiltonianCycle.push_back(vertex);
            prevVertex = vertex;
        }
    }

    // Adding the edge to complete the cycle
    if (!hamiltonianCycle.empty() && prevVertex != -1 && hamiltonianCycle.front() != prevVertex) {
        int firstVertex = hamiltonianCycle.front();
        for (const auto& [w, weight] : graph.getAdjacencyList()[prevVertex]) {
            if (w == firstVertex) {
                totalWeight += weight;
                break;
            }
        }
    }

    return {hamiltonianCycle, totalWeight};
}

// Prim's algorithm implementation to compute the MST of a graph
MSTResult GraphAlgorithm::primMST(const Graph& graph) {
    int V = graph.getNumberOfVertices(); // Get the number of vertices in the graph
    std::vector<bool> inMST(V, false); // Vector to keep track of vertices in the MST, initially all false
    std::vector<double> key(V, std::numeric_limits<double>::max()); // Key values to pick the minimum weight edge
    std::vector<int> parent(V, -1); // Array to store the constructed MST

    // Priority queue to pick the minimum weight edge at every step
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;

    int startVertex = 0;
    pq.push({0.0, startVertex});
    key[startVertex] = 0.0;

    MSTResult mstResult;

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        if (inMST[u]) continue;
        inMST[u] = true;

        if (parent[u] != -1) {
            mstResult.addMSTEdge(parent[u], u, key[u]);
        }

        const auto& adj = graph.getAdjacencyList();
        for (const auto& [v, weight] : adj[u]) {
            if (!inMST[v] && weight < key[v]) {
                key[v] = weight;
                pq.push({weight, v});
                parent[v] = u;
            }
        }
    }

    return mstResult;
}

// Kruskal's algorithm implementation to compute the MST of a graph
MSTResult GraphAlgorithm::kruskalMST(const Graph& graph) {
    std::vector<std::tuple<double, int, int>> edges; // weight, start vertex, end vertex
    const auto& adjList = graph.getAdjacencyList(); // Get the adjacency list of the graph

    for (int v = 0; v < adjList.size(); ++v) { 
        for (const auto& edge : adjList[v]) {
            if (v < edge.first) { // Only consider edges once (undirected graph), ensure each edge is only added once
                edges.emplace_back(edge.second, v, edge.first);
            }
        }
    }

    // Sort edges by ascending weight
    std::sort(edges.begin(), edges.end());

    UnionFind uf(graph.getNumberOfVertices()); // Initialize a union-find data structure with the number of vertices
    MSTResult mstResult;

    for (const auto& [weight, u, v] : edges) {
        if (uf.find(u) != uf.find(v)) { // Check if adding the edge creates a cycle
            uf.unite(u, v); // If not, add the edge to the MST
            mstResult.addMSTEdge(u, v, weight); // Add the edge to the MST result
            if (mstResult.getEdges().size() == graph.getNumberOfVertices() - 1) break; // Stop when MST is complete
        }
    }

    return mstResult;
}

// Count the number of connected components in the graph using DFS
int GraphAlgorithm::countConnectedComponents(const Graph& graph) {
    int V = graph.getNumberOfVertices();
    std::vector<bool> visited(V, false);
    int count = 0;

    for (int i = 0; i < V; i++) {
        if (!visited[i]) {
            dfs(i, visited, graph);
            count++;
        }
    }

    return count;
}

// Depth-first search to explore the graph
void GraphAlgorithm::dfs(int v, std::vector<bool>& visited, const Graph& graph) {
    visited[v] = true;
    for (const auto& edge : graph.getAdjacencyList()[v]) {
        if (!visited[edge.first]) {
            dfs(edge.first, visited, graph);
        }
    }
}
