#include "GraphAlgorithm.h"
#include <queue>
#include <limits>
#include <iostream>
#include <stdexcept>

std::pair<std::vector<double>, std::vector<int>> GraphAlgorithm::dijkstra(const Graph& graph, int source) 
{
    int n = graph.numVertices();
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<int> pred(n, -1); // Predecessor array
    dist[source] = 0;

    using P = std::pair<double, int>; 
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq; 
    pq.push({0, source});

    while (!pq.empty()) 
    {
        double d = pq.top().first; // Current distance
        int u = pq.top().second; // Current vertex
        pq.pop(); // Remove

        if (d > dist[u]) 
        {
            continue;
        }

        // Relaxation step: for each edge (u, v), if dist[u] + weight < dist[v], update dist[v]
        for (const auto& edge : graph.getEdges(u)) 
        {
            int v = edge.target;
            double weight = edge.weight;
            if (dist[u] + weight < dist[v]) 
            {
                dist[v] = dist[u] + weight; 
                pred[v] = u;
                pq.push({dist[v], v});
            }
        }
    }

    return {dist, pred};
}

std::pair<std::vector<double>, std::vector<int>> GraphAlgorithm::mooreBellmanFord(const Graph& graph, int source) 
{
    int n = graph.numVertices();
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<int> pred(n, -1);
    dist[source] = 0;

    // Relaxing edges n - 1 times --> it ensures that the shortest path is found, if there is no negative weight cycle
    for (int i = 1; i < n; ++i) 
    {
        for (int u = 0; u < n; ++u) 
        {
            for (const auto& edge : graph.getEdges(u)) 
            {
                int v = edge.target;
                double weight = edge.weight;
                if (dist[u] + weight < dist[v]) 
                {
                    dist[v] = dist[u] + weight;
                    pred[v] = u;
                }
            }
        }
    }

    // Checking for negative weight cycles, if any edge can be relaxed, then a negative weight cycle exists
    for (int u = 0; u < n; ++u) 
    {
        for (const auto& edge : graph.getEdges(u)) 
        {
            int v = edge.target;
            double weight = edge.weight;
            if (dist[u] + weight < dist[v]) 
            {
                throw std::runtime_error("Graph contains a negative weight cycle");
            }
        }
    }

    return {dist, pred};
}
