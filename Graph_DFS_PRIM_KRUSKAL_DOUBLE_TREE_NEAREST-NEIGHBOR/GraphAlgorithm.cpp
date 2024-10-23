#include "GraphAlgorithm.h"
#include <limits>
#include <vector>
#include <queue>
#include <set>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <stdexcept> 

// Nearest Neighbor TSP with total weight
std::pair<std::vector<int>, double> GraphAlgorithm::nearestNeighborTSP(const Graph& graph) 
{
    int n = graph.getNumberOfVertices(); 
    std::vector<bool> visited(n, false); // Initialize all vertices as not visited
    std::vector<int> tour; // Store the TSP tour
    int current = 0; // Start from vertex 0
    double totalWeight = 0.0;

    tour.push_back(current); // Add the starting vertex to the tour
    visited[current] = true; // Mark the starting vertex as visited

    for (int i = 1; i < n; ++i) // Visit the remaining n-1 vertices, because we already visited the first vertex
    {
        /*
        * On the first comparison, any distance will be less than the maximum value and min_distance will be updated.
        * On subsequent comparisons, if the distance is less than the current min_distance, it will be updated.
        */
        double min_distance = std::numeric_limits<double>::max(); //provides the maximum value for the data type
        int next = -1;

        for (int j = 0; j < n; ++j) // Finds the nearest unvisited neighbor of the current vertex
        {
            if (!visited[j]) // Check if the vertex has not been visited
            {
                try {
                    double distance = graph.getEdgeWeight(current, j); // Get the weight of the edge between the current vertex and the neighbor
                    if (distance < min_distance) 
                    {
                        min_distance = distance; // Update the minimum distance
                        next = j; // Update the next vertex to visit
                    }
                } catch (const std::out_of_range& e) 
                {
                    // Edge does not exist, skip
                }
            }
        }

        totalWeight += min_distance; // Add the weight of the edge to the total weight
        current = next; // Move to the next vertex, j becomes the current vertex
        tour.push_back(current); // Add the next vertex to the tour
        visited[current] = true; // Mark the next vertex as visited
    }

    // Close the tour by returning to the starting node
    totalWeight += graph.getEdgeWeight(current, tour[0]); 
    tour.push_back(tour[0]);

    return {tour, totalWeight};
}

// Helper function to create a doubled graph (Double-Tree)
Graph GraphAlgorithm::doubleEdges(const Graph& graph) 
{
    Graph doubledGraph(graph.getNumberOfVertices()); // Create a new graph with the same number of vertices
    const auto& adj = graph.getAdjacencyList(); // Get the adjacency list of the original graph

    for (int v = 0; v < graph.getNumberOfVertices(); ++v) 
    {
        for (const auto& [w, weight] : adj[v]) // Iterate over the adjacency list of each vertex
        {
            doubledGraph.addEdge(v, w, weight); 
            doubledGraph.addEdge(w, v, weight);  // Add both directions
        }
    }

    return doubledGraph;
}

// Helper function to find an Eulerian tour (Double-Tree)
std::vector<int> GraphAlgorithm::findEulerianTour(const Graph& graph) 
{
    /*
    * Use a stack to traverse the graph
    * Push the starting vertex onto the stack
    * While the stack is not empty:
    * - If the current vertex has no neighbors, add it to the tour and pop it from the stack
    * - Otherwise, push an adjacent vertex onto the stack and remove the edge between them
    */
    std::vector<int> tour; // Store the Eulerian tour
    std::stack<int> stack; // Stack to keep track of the vertices, for DFS traversal
    /*
    * Unordered map where each key is a vertex and the value is a multiset of adjacent vertices.
    * The miltiset in 'adjList' allows for handling multiple edges between the same pair of vertices if present.
    * This list will be used to keep track of the edges in the graph and remove them as we traverse.
    */
    std::unordered_map<int, std::multiset<int>> adjList; 

    const auto& adj = graph.getAdjacencyList(); //
    for (int v = 0; v < graph.getNumberOfVertices(); ++v) 
    {
        for (const auto& [w, weight] : adj[v]) 
        {
            adjList[v].insert(w); // Populate the 'adjList' with all edges from the original graph
        }
    }

    /* DSF Traversal */

    int start = 0;  // Starting vertex
    stack.push(start); // Push the starting vertex onto the stack

    while (!stack.empty()) 
    {
        int v = stack.top(); // Get the current vertex from the top of the stack
        if (adjList[v].empty()) // Check if v has any remaining edges in adjacency list
        {
            tour.push_back(v); // v is added to the tour, since there are no more edges to explore from it
            stack.pop(); // The stack is popped to backtrack to the previous vertex
        } else 
        {
            /* 
            * If v has remaining edges:
            * - Get the first adjacent vertex u from the multiset of adjacent vertices of v
            *   adjsList[v].begin() returns an iterator to the first element in the multiset of adjacent vertices of v, which is u.
            */
            int u = *adjList[v].begin(); 
            stack.push(u);// Push u onto the stack

            // Remove the edge between v and u from the adjacency list, this ensures that the edge is visited only once
            adjList[v].erase(adjList[v].begin()); 
            adjList[u].erase(adjList[u].find(v)); 
        }
    }

    return tour;
}

// Double Tree Algorithm for TSP approximation
std::pair<std::vector<int>, double> GraphAlgorithm::doubleTreeTSP(const Graph& graph) 
{
    // Step 1: Construct the MST using Prim's or Kruskal's algorithm
    MSTResult mstResult = kruskalMST(graph);
    Graph mstGraph(graph.getNumberOfVertices());
    for (const auto& edge : mstResult.edges) // Iterate over the edges of the MST, each edge is a pair of vertices (u, v)
    {
        int u = edge.first; 
        int v = edge.second;
        // Find the weight of the edge (u, v) in the original graph
        double weight = -1; // To ensure that the weight is found!
        for (const auto& [w, edgeWeight] : graph.getAdjacencyList()[u]) // iterate over the adjacency list of vertex u in the original graph, returning a list of neighbors of u and their weights
        {
            if (w == v) // For each neighbor w of u, with an edge weight edgeWeight, check if w is equal to v, meaning that the edge (u, v) exists in the original graph
            {
                weight = edgeWeight; // The weight of the edge (u, v) is found and we assign it to the variable weight
                break;
            }
        }
        mstGraph.addEdge(u, v, weight); // Add the edge (u, v) with weight to the MST graph
    }

    // Step 2: Double the edges of the MST
    Graph doubledMST = doubleEdges(mstGraph);

    // Step 3: Compute an Eulerian tour of the doubled MST
    std::vector<int> eulerianTour = findEulerianTour(doubledMST);

    // Step 4: Transform the Eulerian tour into a Hamiltonian cycle
    std::vector<int> hamiltonianCycle; // Store the vertices of the Hamiltonian cycle
    std::unordered_set<int> visited; // Keep track of visited vertices, each element is unique in the set
    double totalWeight = 0.0;

    int prevVertex = -1; // Keep track of the previous vertex, initialized to -1 to indicate the starting vertex
    for (int vertex : eulerianTour) // Iterate over the vertices in the Eulerian tour
    {
        if (visited.find(vertex) == visited.end()) // Checks if the current vertex has already been visited
        {
            if (prevVertex != -1) // Ensures that there is a valid previous vertex before attempting to find and add the edge weight between prevVertex and the current vertex
            {
                // Find the weight of the edge (prevVertex, vertex) in the original graph
                for (const auto& [w, weight] : graph.getAdjacencyList()[prevVertex]) 
                {
                    if (w == vertex) 
                    {
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
    if (!hamiltonianCycle.empty() && prevVertex != -1 && hamiltonianCycle.front() != prevVertex) // Check if the Hamiltonian cycle is not empty and the last vertex is not the same as the first vertex
    {
        int firstVertex = hamiltonianCycle.front(); // Get the first vertex of the Hamiltonian cycle
        for (const auto& [w, weight] : graph.getAdjacencyList()[prevVertex]) // Iterate over the adjacency list of the previous vertex to finf the edge connecting the last vertex to the first vertex in the Hamiltonian cycle
        {
            if (w == firstVertex) 
            {
                totalWeight += weight;
                break;
            }
        }
    }

    return {hamiltonianCycle, totalWeight};
}

// Prim's algorithm implementation to compute the MST of a graph
MSTResult GraphAlgorithm::primMST(const Graph& graph) 
{
    int V = graph.getNumberOfVertices();
    std::vector<bool> inMST(V, false); //Keep track of vertices in the MST, initialized to false
    std::vector<double> key(V, std::numeric_limits<double>::max()); //Vector to store the minimum weight edge that connects each vertex to the MST, initialized to the maximum value
    std::vector<int> parent(V, -1); //Vector to store the parent of each vertex in the MST, initialized to -1

    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq; //Min-heap to store the vertices with the minimum key value, ordered by key value

    int startVertex = 0; // Start from vertex 0
    pq.push({0.0, startVertex}); // Push the starting vertex onto the pq with key value 0
    key[startVertex] = 0.0;

    MSTResult mstResult; //store the edges of the MST

    /*Chek the Priority Queue*/
    while (!pq.empty()) 
    {
        /*Extract the minimum key vertex*/
        int u = pq.top().second; // Get the vertex with the minimum key value .second (second element of the pair, which is the vertex itself, the first element is the key value)
        pq.pop(); // Remove the vertex from the pq

        /*Check inclusion in MST*/
        if (inMST[u]) continue; // Skip if the vertex is already in the MST
        inMST[u] = true; // Mark the vertex as part of the MST

        /*Add edge to mstResult*/
        if (parent[u] != -1) // If u has a valid parent
        {
            mstResult.addMSTEdge(parent[u], u, key[u]); //add the edge from the parent to u to the MST with the weight key[u]
        }

        /*Update Adjacent Vertices*/
        const auto& adj = graph.getAdjacencyList(); 
        for (const auto& [v, weight] : adj[u]) // Get the adj list of vertex u. For each adj vertex v and the corresponding edge weight
        {
            if (!inMST[v] && weight < key[v]) // v is not in the MST and the weight of the edge (u, v) is less than the current key value of v
            {
                key[v] = weight; // Update the key value of v to the weight of the edge (u, v)
                pq.push({weight, v}); // Push the updated key value of v and v onto the pq
                parent[v] = u; // Update the parent of v to u
            }
        }
    }

    return mstResult;
}

// Kruskal's algorithm implementation to compute the MST of a graph
MSTResult GraphAlgorithm::kruskalMST(const Graph& graph) 
{
    /*Initialize the edge list and get adjList*/
    std::vector<std::tuple<double, int, int>> edges; // Vector to store all edges of the graph in the form of tuples (weight, u, v).
    const auto& adjList = graph.getAdjacencyList(); // Reference to the adjacency list of the graph

    /*Populate the edge list*/
    for (int v = 0; v < adjList.size(); ++v) // Iterate through each vertex v in the adjacency list
    {
        for (const auto& edge : adjList[v]) // For each vertex v, iterate through its neighbors
        {
            if (v < edge.first) // Add the edge only once, either (u, v) or (v, u) but not both, an edge (u,v) appears twice in the adj list of u and v
            {
                edges.emplace_back(edge.second, v, edge.first); //Add the edge to the 'edges' vector as a tuple (weight, u, v)
            }
        }
    }

    /*Sort the edges by weight*/
    std::sort(edges.begin(), edges.end()); // Default tuple comparison is based on the first element, which is the weight, then u, then v

    /*Initialize the Union-Find Structure and MSTResul*/
    UnionFind uf(graph.getNumberOfVertices()); // Create a Union-Find data structure with the number of vertices in the graph
    MSTResult mstResult;

    /*Process the edges*/
    for (const auto& [weight, u, v] : edges) // Iterate over each edge (weight, u, v) in the sorted edges vector
    {
        if (uf.find(u) != uf.find(v)) // Check for cycles, check if u and v are in the same connected component
        {
            uf.unite(u, v); // If u and v are not in the same connected component, unite them
            mstResult.addMSTEdge(u, v, weight); // Add it to the MST result
            if (mstResult.getEdges().size() == graph.getNumberOfVertices() - 1) break; // Stop when the MST has V-1 edges, otherwise it will be a cycle
        }
    }

    return mstResult;
}

// Count the number of connected components in the graph using DFS
int GraphAlgorithm::countConnectedComponents(const Graph& graph) 
{
    /*Initialize the visited vector and count*/
    int V = graph.getNumberOfVertices();
    std::vector<bool> visited(V, false);
    int count = 0;

    /*DFS to explore the graph*/
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
    /*Mark the vertex as visited and explore its neighbors*/
    visited[v] = true;
    /*Iterate over the neighbors of the vertex*/
    for (const auto& edge : graph.getAdjacencyList()[v]) 
    {
        if (!visited[edge.first]) 
        {
            dfs(edge.first, visited, graph);
        }
    }
}