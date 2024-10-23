#include "algo.h"

std::vector<bool> depth_first_search(Graph& graph, int start, std::vector<bool>& visited) {
    std::stack<int> stack;
    stack.push(start);
    
    visited[start] = true;

    while (!stack.empty()) {
        int current = stack.top();
        stack.pop();
        // std::cout << current << " ";

        for (std::pair<int, double> neighbor : graph.get_neighbor(current)) {
            if (!visited[neighbor.first]) {
                stack.push(neighbor.first);
                visited[neighbor.first] = true;
            }
        }
    }
    return visited;
}

std::vector<int> breadth_first_search_Max_Flow(Graph& graph, int start, std::vector<bool>& visited) {
    std::list<int> list;
    std::vector<int> order;
    list.push_back(start);
    
    visited[start] = true;
    order.push_back(start);

    while (!list.empty()) {
        int current = list.front();
        list.pop_front();
        // std::cout << current << " ";

        for (std::pair<Node, kostMinEdge>& neighbor : graph.get_neighbor_KostMin(current)) {
            if (!visited[neighbor.second.node2.node1]) {
                list.push_back(neighbor.second.node2.node1);
                visited[neighbor.second.node2.node1] = true;
                order.push_back(neighbor.second.node2.node1);
            }
        }
    }

    return order;
}

std::vector<std::vector<int> > component_bfs(Graph& graph) {
    std::vector<std::vector<int> > orders;
    std::vector<bool> visited(graph.get_V(), false);

    for (int i = 0; i < graph.get_V();i++) {
        if (!visited[i]) {
            orders.push_back(breadth_first_search_Max_Flow(graph, i, visited));
        }
    }
    return orders;
}

int connected_components(Graph& graph) {
    int components = 0;
    std::vector<bool> visited(graph.get_V(), false);

    for (int i = 0; i < graph.get_V(); i++) {
        if (visited[i] == false) {
            visited = depth_first_search(graph, i, visited);
            components++;
        }
    }
    return components;
}

std::vector<int> breadth_first_search_to_target(Graph& graph, int source, int target, double& kleinsterFluss) {
    int V = graph.get_V();
    std::vector<bool> visited(V, false);
    std::list<int> list;
    std::vector<int> order;
    std::vector<int> parent(V, -1);
    double kf = 0.0;


    list.push_back(source);
    visited[source] = true;
    parent[source] = source;

    while (!list.empty() && !visited[target]) {
        int current = list.front();
        list.pop_front();

        for (std::pair<int, double>& neighbor : graph.get_neighbor(current)) {
            int node = neighbor.first;
            double weight = neighbor.second;
            if (!visited[node] && weight != 0) {
                parent[node] = current;
                list.push_back(node);
                visited[node] = true;
            }
        }
    }
    if (visited[target]) {
        order.push_back(target);
        int current = parent[target];
        kf = graph.get_Weight(current, target);
        while (current != source) {
            order.push_back(current);
            kf = std::min(kf, graph.get_Weight(parent[current], current));
            current = parent[current];
        }
        order.push_back(source);
        kleinsterFluss = kf;
    }
    return order;
}

// Breadth-First-Search for Cost-Minimal Flow Graphs
std::vector<int> breadth_first_search_to_target_KostMinFlow(Graph& graph, int source, int target, double& kleinsterFluss) {
    int V = graph.get_V() + 2;
    std::vector<bool> visited(V, false);
    std::list<int> list;
    std::vector<int> order;
    std::vector<int> parent(V, -1);
    double kf = 0.0;

    list.push_back(source);
    visited[source] = true;
    parent[source] = source;

    while (!list.empty() && !visited[target]) {
        int current = list.front();
        list.pop_front();

        for (std::pair<Node, kostMinEdge>& neighbor : graph.get_neighbor_KostMin(current)) {
            int node = neighbor.second.node2.node1;
            double weight = neighbor.second.kapazitaet;
            if (!visited[node] && weight != 0) {
                parent[node] = current;
                list.push_back(node);
                visited[node] = true;
            }
        }
    }
    if (visited[target]) {
        order.push_back(target);
        int current = parent[target];
        kf = graph.get_Kapa(current, target);
        while (current != source) {
            order.push_back(current);
            kf = std::min(kf, graph.get_Kapa(parent[current], current));
            current = parent[current];
        }
        order.push_back(source);
        kleinsterFluss = kf;
    }
    return order;
}

// Edmonds-Karp Algorithm for Cost-Minimal Flow Graphs
double edmonds_karp_KostMinFlowGraph(Graph& graph, int source, int target) {
    
    graph.addResidualEdgesMaxFlow();
    
    double maximalerFluss = 0.0;
    std::vector<int> order;
    double kleinsterFluss = std::numeric_limits<double>::infinity();

    while (true) {
        order = breadth_first_search_to_target_KostMinFlow(graph, source, target, kleinsterFluss);
        if (order.empty()) {
            break;
        } else {
            for (int i = order.size()-1; i > 0; i--) {
                graph.adjust_Kapa(order[i], order[i-1], -kleinsterFluss);
                graph.adjust_Fluss(order[i], order[i-1], kleinsterFluss);
                graph.adjust_Kapa(order[i-1], order[i], kleinsterFluss);
                graph.adjust_Fluss(order[i-1], order[i], kleinsterFluss);
            }
            maximalerFluss += kleinsterFluss;
        }
    }
    
    return maximalerFluss;
}

// Add Super Source and Sink 
// Idea: To transform the proble into a max flow problem
void addSuperSourceAndSink(Graph& graph) {
    int V = graph.get_V();
    Node sourceNode;
    Node sinkNode;
    sourceNode.node1 = V;
    sinkNode.node1 = V+1;
    sourceNode.balance = 0;
    sinkNode.balance = 0;

    for (Node& node : graph.get_NodeList()) {
        if (node.balance > 0) {
            kostMinEdge edge = {node, 0, 0,node.balance};
            graph.addEdgeToListKostMin(sourceNode.node1, sourceNode, edge);
        } else if (node.balance < 0) {
            kostMinEdge edge = {sinkNode, 0, 0,-node.balance};
            graph.addEdgeToListKostMin(node.node1, node, edge);
        } else {}
    }
}

void removeSuperSourceAndSink(Graph& graph) {
    graph.removeNodeKostMin();
}

// Moore-Bellman-Ford Algorithm for Negative Cycles
std::vector<int> moore_bellman_ford_negative_cycle(Graph& graph) {
    double INF = std::numeric_limits<double>::infinity();
    int V = graph.get_V();

    std::vector<double> dist;  
    std::vector<int> prev;  
    std::vector<int> cycle;

    // Breath-First-Search
    std::vector<std::vector<int> > orders = component_bfs(graph);

    for (std::vector<int>& order : orders) {
        
        dist.assign(V, INF);
        prev.assign(V, -1);
        cycle.clear();

        dist[order[0]] = 0;  
        prev[order[0]] = order[0];

        for (int i = 0; i < order.size()-1; i++) {
            // Over all edges iteration
            for (int& j : order) {
                for (auto neighbour : graph.get_neighbor_KostMin(j)){
                    int node1 = j;
                    int node2 = neighbour.second.node2.node1;
                    double kosten = neighbour.second.kosten;

                    if (dist[node1] != INF && dist[node1] + kosten < dist[node2]) {
                        dist[node2] = dist[node1] + kosten;
                        prev[node2] = node1;
                    }
                }
            }
        }

        // Negative Cycle Detection
        for (int& j : order) {
            for (auto& neighbour : graph.get_neighbor_KostMin(j)){
                int node1 = j;
                int node2 = neighbour.second.node2.node1;
                double kosten = neighbour.second.kosten;
                if (dist[node1] != INF && dist[node1] + kosten < dist[node2]) {
                    int current = node2;
                    for (int n = 0; n < order.size(); n++) {
                        current = prev[current];
                    }
                    int loopNode = current;
                    while (true) {
                        cycle.push_back(current);
                        current = prev[current];
                        if (current == loopNode)
                            break;
                    }
                    cycle.push_back(loopNode);
                    std::reverse(cycle.begin(), cycle.end());
                    return cycle;
                }
            }
        }
    }
    return cycle;
}

// Moore-Bellman-Ford Algorithm for Cost-Minimal Flow Graphs
std::vector<int> mbf(Graph& graph, int& source, int& sink) {
    double INF = std::numeric_limits<double>::infinity();
    int V = graph.get_V();

    std::vector<bool> visited(V, false);
    std::vector<double> dist(V, INF);  
    std::vector<int> prev(V, -1);  
    std::vector<int> path;

    // Breath-First-Search
    std::vector<int> order = breadth_first_search_Max_Flow(graph, source, visited);


    dist[source] = 0;  
    prev[source] = source;

    for (int i = 0; i < order.size()-1; i++) {
        // Over all edges iteration
        for (int& j : order) {
            for (auto neighbour : graph.get_neighbor_KostMin(j)){
                int node1 = j;
                int node2 = neighbour.second.node2.node1;
                double kosten = neighbour.second.kosten;

                if (dist[node1] != INF && dist[node1] + kosten < dist[node2]) {
                    dist[node2] = dist[node1] + kosten;
                    prev[node2] = node1;
                }
            }
        }
    }

    // Negative Cycle Detection
    for (int& j : order) {
        for (auto& neighbour : graph.get_neighbor_KostMin(j)){
            int node1 = j;
            int node2 = neighbour.second.node2.node1;
            double kosten = neighbour.second.kosten;
            if (dist[node1] != INF && dist[node1] + kosten < dist[node2]) {
                std::cout << "Negative Cycle found!" << std::endl;
                return path;
            }
        }
    }

    if (prev[sink] != -1) {
        int current = sink;
        while (current != source) {
            path.push_back(current);
            current = prev[current];
        }
        path.push_back(source);
        std::reverse(path.begin(), path.end());
        return path;
    }

    return path;

}

//
bool graph_Loesbarkeit(Graph& graph, double maxFluss) {
    // Suppy and Demand calculation --> supply = demand & demand = nachfrage
    double angebot = 0.0;
    double nachfrage = 0.0;
    for (Node& node : graph.get_NodeList()) {
        if (node.balance < 0) {
            nachfrage += node.balance;
        } else if (node.balance > 0) {
            angebot += node.balance;
        } else {}
    }
    //Debbugging
    //std::cout << "Angebot: " << angebot << std::endl;
    //std::cout << "Nachfrage: " << -nachfrage << std::endl;

    // Verify if it is possible to have a b-Fluss
    if(maxFluss != -nachfrage || angebot != -nachfrage) {
        std::cout << "Graph nicht lösbar! Kein b-Fluss Möglich!" << std::endl;
        return false;
    } else {
        std::cout << "Graph lösbar! b-Fluss gefunden!" << std::endl;
        return true;
    }
}

// Update Flow
void updateFlow(Graph& originalGraph, Graph& graph) {
    for (auto& adj : graph.get_adjListKostMin()) {
        for(auto& edge : adj) {
            originalGraph.adjust_Fluss(edge.first.node1, edge.second.node2.node1, edge.second.fluss);
        }
    }
}

// Cycle Canceling Algorithm *********************************
double cycle_canceling(Graph& graph) {
    // Input: Gerichteter Graph G = (V, E) mit oberen Kapazitäten u, Kantenkosten c, Balance b
    // Output: Fluss f mit minimalen Kosten

    // Schritt 1: Überprüfen der Lösbarkeit und Berechnen eines b-Flusses
    Graph operationalGraph = graph; // Copy of the original graph
    Graph maxFlussGraph = graph; // Copy of the original graph

    addSuperSourceAndSink(maxFlussGraph);

    // MaxFluss calculation using Edmonds-Karp Algorithm
    double maxFluss = edmonds_karp_KostMinFlowGraph(maxFlussGraph, maxFlussGraph.get_adjListKostMin().size()-2, maxFlussGraph.get_adjListKostMin().size()-1);
    
    // Remove Super Source and Sink to restore the original graph
    removeSuperSourceAndSink(maxFlussGraph);

    // Feasilibility check
    bool loesbar = graph_Loesbarkeit(maxFlussGraph, maxFluss);
    //If the max flow is equal to the total supply/demand) then transfer flow to original graph, update flow
    if (loesbar) {
        updateFlow(operationalGraph, maxFlussGraph);
    } else {
        return 0.0;
    }

    // Schritt 2: Bestimmen des ResidualGraphen
    Graph residualGraph;
    std::vector<int> cycle;
    double kleinsterFluss;
    double kostenMinimalerFluss = 0.0;

    // Schritt 2: Bestimmung des Residualgraphen
    while (true) {
        
        residualGraph = operationalGraph; 
        residualGraph.addResidualEdgesCostKapa();

        // Schritt 3: Negativen Zykel with Moore-Bellman-Ford Algorithmus 
        cycle = moore_bellman_ford_negative_cycle(residualGraph);
        

        if (cycle.empty()) {
            std::cout << "No negative cycle found!" << std::endl;
            break;
        }

        // Schritt 4: Flussaugmentierung entlang des negativen Zykels
        kleinsterFluss = std::numeric_limits<double>::infinity(); // Set to infinity to detect any smaller capacity will be found

        for (int i = 1; i < cycle.size(); i++) { // Find the smallest capacity in the cycle
            kleinsterFluss = std::min(kleinsterFluss, residualGraph.get_Kapa(cycle[i-1], cycle[i])); // Retrieves the capacity of the edge in the residual graph
        }
        
        for (int i = 1; i < cycle.size(); i++) {
            operationalGraph.adjust_Fluss(cycle[i-1], cycle[i], kleinsterFluss); // Increase the flow along the edge
            operationalGraph.adjust_Fluss(cycle[i], cycle[i-1], -kleinsterFluss); // Decrease the flow along the edge
        }
        kostenMinimalerFluss = operationalGraph.get_GesamtKosten(); // Update the total cost of the flow, getGeamtKosten() calculates the total cost of the current flow
    
    }

    return kostenMinimalerFluss; // Total cost o the minimum cost flow
}

// Balancen ausgleichen
bool balancen_ausgelichen(std::vector<Node> nodeList, std::vector<double>& new_balance) {
    // Loops over all nodes and checks if the balance of the node is equal to the new balance
    for (int i = 0; i < nodeList.size(); i++) {
        if (nodeList[i].balance != new_balance[i]) {
            return false;
        }
    }
    return true;
}

// Successive Shortest Path Algorithm **********************************
double successive_shortest_Path(Graph& graph) {

    // Schritt 1: Konstruktion kostenminimaler Fluss mit beliebigen Balance
    // Idee: Kanten mit negativen Kosten werden mit maximalen Flusswert auslasten
    int V = graph.get_V();

    // Maximize flow on edges with negative costs
    // Ensures that edges with negative costs are fully utilized to minimize initial costs.
    for (auto& adj : graph.get_adjListKostMin()) {
        for (auto& edge : adj) {
            if (edge.second.kosten < 0) {
                graph.adjust_Fluss(edge.first.node1, edge.second.node2.node1, edge.second.kapazitaet);
            }
        }
    }

    // Schritt 2: Annäherung der Balance b_strich an b, b_strich keeps track of the net flow into and out of each node
    
    std::vector<double> b_strich(V, 0.0); // Initialize the balance vector
    for (auto& adj : graph.get_adjListKostMin()) { // iterate over the adj list
        for (auto& edge : adj) { // iterate over the edges
            b_strich[edge.first.node1] += edge.second.fluss; // Increase the balance of the starting node of the edge
            b_strich[edge.second.node2.node1] -= edge.second.fluss; // Decrease the balance of the ending node of the edge
        }
    }

    std::vector<int> sources;
    std::vector<int> sinks;
    Graph residualGraph;
    std::vector<int> weg;
    std::vector<int> kuerzesterWeg; // Shortest Path Vector
    double kleinsterFluss;
    int tmp_Source = 0;
    int tmp_Sink = 0;


    // Iterates until the balance is balanced
    while (!balancen_ausgelichen(graph.get_NodeList(), b_strich)){
        weg.clear();
        kuerzesterWeg.clear();
        kleinsterFluss = std::numeric_limits<double>::infinity(); // Max flow that can be pushed along the path
        sources.clear();
        sinks.clear();
        
        // Identify sources and sinks
        for (int i = 0; i < V; i++) {
            if ((graph.get_NodeList()[i].balance - b_strich[i]) > 0) {
                sources.push_back(i);
            } else if ((graph.get_NodeList()[i].balance - b_strich[i]) < 0) {
                sinks.push_back(i);
            } else {}
        }

        // Construct the residual graph
        residualGraph = graph;
        residualGraph.addResidualEdgesCostKapa();

        // Path from Source to Sink
        for (int source : sources) {
            for (int sink : sinks) {
                weg = breadth_first_search_to_target_KostMinFlow(residualGraph, source, sink, kleinsterFluss);
                if (!weg.empty()) {
                    tmp_Source = source;
                    tmp_Sink = sink;
                    break;
                }
            }
        }

        if (weg.empty()) {
            std::cout << "Kein Weg gefunden!" << std::endl;
            return 0.0;
        }
        
        // Schritt 3: Suche des kürzesten Weges im Residualgraphen
        // Idee: Moore-Bellman-Ford Algorithmus to find shortest path with minimum cost in residual graph
        
        kuerzesterWeg = mbf(residualGraph, tmp_Source, tmp_Sink); // Find Shortest Path with Moore-Bellman-Ford Algorithm


        // Schritt 4: Bestimmung von y als Minimum der Werte: min(Minimale Residualkapazität der Kanten aus p, Maximaler Wert den b_strich(s) steigen darf, Maximaler Wert den b_strich(t) sinken darf)
        for (int i = 1; i < kuerzesterWeg.size(); i++) {
            kleinsterFluss = std::min(kleinsterFluss, residualGraph.get_Kapa(kuerzesterWeg[i-1], kuerzesterWeg[i])); // Find the smallest capacity in the path
        }
       
       // Schritt 5: b-Flussanpassung entlang des Weges p um Wert y
        kleinsterFluss = std::min(kleinsterFluss, graph.get_NodeList()[tmp_Source].balance - b_strich[tmp_Source]); // Limit flow by source balance (difference source node original balance and its current imbalance)
        
        kleinsterFluss = std::min(kleinsterFluss, b_strich[tmp_Sink] - graph.get_NodeList()[tmp_Sink].balance); // Limit flow by sink balance (dif. current imbalance and sink node original balance)
        
        // Adjust flow along the path
        for (int i = 1; i < kuerzesterWeg.size(); i++) {
            graph.adjust_Fluss(kuerzesterWeg[i-1], kuerzesterWeg[i], kleinsterFluss); // Increases the flow along the edge, foward direction
            graph.adjust_Fluss(kuerzesterWeg[i], kuerzesterWeg[i-1], -kleinsterFluss); // Decreases the flow along the edge, backward direction

            // Balance b_strich anpassen um den Wert y
            b_strich[kuerzesterWeg[i-1]] += kleinsterFluss; // Increase the balance 
            b_strich[kuerzesterWeg[i]] -= kleinsterFluss; // Decrease the balance
        }
    }
    return graph.get_GesamtKosten();
}