#include "GraphAlgorithm.h"
#include <iostream>

/**
 * Calculates the cost of a given tour in a graph.
 * 
 * @param graph The graph containing the edges and weights.
 * @param tour The tour is a vector representing the sequence of vertices visited
 * @return The cost of the tour.
 */
double GraphAlgorithm::calculateTourCost(const Graph& graph, const std::vector<int>& tour) 
{
    double cost = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i) // Iterate over the tour and sum the weights of the edges between consecutive vertices
    {
        cost += graph.getEdgeWeight(tour[i], tour[i + 1]); 
    }
    cost += graph.getEdgeWeight(tour.back(), tour.front()); // Add the weight of the last vertex back to the first to complete the tour 
    return cost;
}

/**
 * Calculates the minimum cost of visiting all vertices in the graph using all possible tours.
 * 
 * @param graph The graph to calculate the tours on.
 * @return The minimum cost of visiting all vertices.
 */
double GraphAlgorithm::allTours(const Graph& graph) 
{
    int V = graph.getNumberOfVertices();
    std::vector<int> vertices(V); 
    std::iota(vertices.begin(), vertices.end(), 0); // Fill the vector with the values from 0 to V-1

    double minCost = std::numeric_limits<double>::infinity(); // Initialize the minimum cost with the maximum value of double

    do 
    {
        double currentCost = calculateTourCost(graph, vertices);
        if (currentCost < minCost) 
        {
            minCost = currentCost;
        }
    } while (std::next_permutation(vertices.begin(), vertices.end())); // Iterate over all permutations of the vertex list

    return minCost;
}

/**
 * Calculates a lower bound on the cost of completing a tour from the current path
 * 
 * @param graph The graph to calculate the lower bound for.
 * @param path The path is a vector of int representing the current partial tour
 * @return The calculated lower bound.
 */
double GraphAlgorithm::calculateLowerBound(const Graph& graph, const std::vector<int>& path) 
{
    int V = graph.getNumberOfVertices();
    std::vector<int> remainingVertices;
    std::vector<bool> inPath(V, false);

    // Identifies the vertices not yet included in the current path
    for (int vertex : path) 
    {
        inPath[vertex] = true;
    }
    for (int i = 0; i < V; ++i) 
    {
        if (!inPath[i]) 
        {
            remainingVertices.push_back(i);
        }
    }

    // If no remaining vertices, return 0, so all vertices are included in the path
    if (remainingVertices.empty()) 
    {
        return 0.0;
    }

    // Finds the minimum weight of an edge connecting the last vertex in the current path to any vertex not in the path
    double minConnectingEdge = std::numeric_limits<double>::infinity(); 
    int lastVertex = path.back();
    for (const auto& [neighbor, weight] : graph.getAdjacencyList()[lastVertex]) 
    {
        if (!inPath[neighbor] && weight < minConnectingEdge) 
        {
            minConnectingEdge = weight;
        }
    }

    return minConnectingEdge; // Returns the minimum weight as the lower bound
}

/**
 * Solves the Traveling Salesman Problem (TSP) using the Branch and Bound algorithm.
 * 
 * @param graph The graph representing the TSP problem.
 * @return The cost of the optimal TSP tour.
 */
double GraphAlgorithm::solveTSPBranchAndBound(const Graph& graph) 
{
    std::priority_queue<Node, std::vector<Node>, CompareNode> pq; // Priority queue to store the nodes in the search
    
    Node root; // Initialize the root node with the first vertex
    root.path.push_back(0); // Start from vertex 0
    root.cost = 0;
    root.lowerBound = calculateLowerBound(graph, root.path); // Calculate the lower bound for the root node
    pq.push(root); // Push the root node into the priority queue

    double bestCost = std::numeric_limits<double>::infinity(); // Initialize the best cost with the maximum value of double

    while (!pq.empty()) 
    {
        Node current = pq.top(); // Get the node with the lowest lower bound from the priority queue
        pq.pop();

        if (current.cost >= bestCost) // Checks if the lower bound of the current node is greater/equal than the best known cost so far
        {
            std::cout << "Pruned node with cost: " << current.cost << " and lower bound: " << current.lowerBound << std::endl;
            continue; // If so, break because the current node cannot lead to a better solution!
        }

        if (current.path.size() == graph.getNumberOfVertices()) // Checks if the current path contains all vertices, i.e., a complete tour
        {
            double totalCost = current.cost + graph.getEdgeWeight(current.path.back(), current.path.front()); // Calculate the total cost of the tour
            if (totalCost < bestCost) // If this total cost is lower than the current best known cost
            {
                bestCost = totalCost; // Update the best known cost
            }
            continue; // Continue to the next node
        }

        // Iterate over all vertices in the graph
        for (int i = 0; i < graph.getNumberOfVertices(); ++i) // Check if the vertex i is already in the current path. if not
        {
            if (std::find(current.path.begin(), current.path.end(), i) == current.path.end()) 
            {
                Node next = current; // A new node is created by copying the current node
                next.path.push_back(i); // The new vertex i is added to the path of next 
                next.cost += graph.getEdgeWeight(current.path.back(), i); // The cost of the edge from the last vertex in the path to i is added to the cost of next
                next.lowerBound = next.cost + calculateLowerBound(graph, next.path); // The lower bound of next is calculated by adding the cost so far to the result of calculateLowerBound
                if (next.lowerBound < bestCost) // If the lower bound of next is less than the best known cost so far, it means that this partial solution is promising, so add it to the pq
                {
                    pq.push(next);
                }
            }
        }
    }

    return bestCost;
}