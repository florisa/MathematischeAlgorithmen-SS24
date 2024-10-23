#include "GraphAlgorithm.h"
#include <queue>
#include <limits>
#include <vector>

bool bfs(const Graph& graph, int source, int sink, std::vector<int>& parent, const std::vector<std::vector<double>>& residualCapacity) {
    std::vector<bool> visited(graph.numVertices(), false);
    std::queue<int> q;
    q.push(source);
    visited[source] = true;
    parent[source] = -1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (const Edge& edge : graph.getEdges(u)) {
            int v = edge.target;
            if (!visited[v] && residualCapacity[u][v] > 0) {
                parent[v] = u;
                visited[v] = true;
                q.push(v);
            }
        }
    }

    return visited[sink];
}

int GraphAlgorithm::edmondsKarp(const Graph& graph, int source, int sink) {
    int n = graph.numVertices();
    std::vector<std::vector<double>> residualCapacity(n, std::vector<double>(n, 0));

    for (int u = 0; u < n; ++u) {
        for (const Edge& edge : graph.getEdges(u)) {
            residualCapacity[u][edge.target] = edge.weight;
        }
    }

    std::vector<int> parent(n);
    int maxFlow = 0;

    while (bfs(graph, source, sink, parent, residualCapacity)) {
        double pathFlow = std::numeric_limits<double>::max();
        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            pathFlow = std::min(pathFlow, residualCapacity[u][v]);
        }

        for (int v = sink; v != source; v = parent[v]) {
            int u = parent[v];
            residualCapacity[u][v] -= pathFlow;
            residualCapacity[v][u] += pathFlow;
        }

        maxFlow += pathFlow;
    }

    return maxFlow;
}
