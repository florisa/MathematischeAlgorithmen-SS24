#include "graph.h"

Graph::Graph() : V(0), adjListKostMin(), adjList(), adjMatrix(), edgeList(), nodeList() {}

Graph::Graph(std::string filename, bool weighted, bool kostMinFluss) : V(0), adjListKostMin(), adjList() ,adjMatrix(), edgeList(), nodeList() {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    if (kostMinFluss) {
        inFile >> V;
        adjListKostMin.resize(V+2);
        int u, v;
        double kosten, kapazitaet;
        double balance;
        std::vector<double> balanceArray(V);
        for(int i = 0; i < V; i++) {
            inFile >> balance;
            Node node = {i, balance};
            nodeList.push_back(node);
            balanceArray[i] = balance;
        }
        while (inFile >> u >> v >> kosten >> kapazitaet) {
            Node node1 = {u, balanceArray[u]};
            Node node2 = {v, balanceArray[v]};
            kostMinEdge edge = {node2, kosten, 0.0,kapazitaet};
            addEdgeToListKostMin(u, node1, edge);
        }
    }
    else if (!kostMinFluss && weighted) {
        inFile >> V;
        int u, v;
        double w;
        while (inFile >> u >> v >> w) {
            addEdge(u, v, w);
        }
    } else {
        inFile >> V;
        int u, v;
        while (inFile >> u >> v) {
            addEdge(u, v, 1);
        }
    }
    inFile.close();
}

Graph::Graph(const Graph& other) : V(other.V), adjListKostMin(other.adjListKostMin), adjList(other.adjList), adjMatrix(other.adjMatrix), edgeList(other.edgeList), nodeList(other.nodeList){}

// Assignment operator
Graph& Graph::operator=(const Graph& other) {
    if (this != &other) {
        V = other.V;
        adjList = other.adjList;
        adjMatrix = other.adjMatrix;
        edgeList = other.edgeList;
        adjListKostMin = other.adjListKostMin;
        nodeList = other.nodeList;
    }
    return *this;
}

void Graph::addEdge(int u, int v, double weight) {
    Edge e(u, v, weight);
    edgeList.push_back(e);
}

void Graph::setAdjList(bool gerichtet) {
    adjList.resize(V);

    if(gerichtet) {
        for (Edge e : edgeList) {
            adjList[e.node1].push_back(std::make_pair(e.node2, e.weight));
        }
    } else {
        for(Edge e : edgeList) {
            adjList[e.node1].push_back(std::make_pair(e.node2, e.weight));
            adjList[e.node2].push_back(std::make_pair(e.node1, e.weight));
        }
    }
}

void Graph::setAdjMatrix(bool gerichtet) {
    adjMatrix.resize(V, std::vector<double>(V, 0.0));

    if(gerichtet) {
        for(Edge e : edgeList) {
            adjMatrix[e.node1][e.node2] = e.weight;
        }

    } else {
        for(Edge e : edgeList) {
            adjMatrix[e.node1][e.node2] = e.weight;
            adjMatrix[e.node2][e.node1] = e.weight;
        }
    }
}

void Graph::addEdgeToList(int u, int v, double weight){
    adjList[u].push_back(std::make_pair(v, weight));
}

void Graph::addEdgeToListKostMin(int u, Node node, kostMinEdge edge) {
    adjListKostMin[u].push_back(std::make_pair(node, edge));
}

void Graph::printAdjList() {
    for (int i = 0; i < V; i++) {
        std::cout << i << ": ";
        for (std::pair<int, double> v : adjList[i]) {
            std::cout << "(" << v.first<< ", " << v.second<< ") ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Graph::printAdjListKostMin() {
    for (auto& i : adjListKostMin) {
        for (std::pair<Node, kostMinEdge>& v : i) {
            std::cout << "(u: " << v.first.node1 << ", u.balance: " << v.first.balance << ", v: " << v.second.node2.node1 << ", v.balance: " << v.second.node2.balance << ", u_v_kosten: " << v.second.kosten << ", u_v_Fluss: " << v.second.fluss << ", u_v_kapazität: " << v.second.kapazitaet << ") ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Graph::printAdjMatrix() {
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            std::cout << adjMatrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int Graph::get_V() {
    return V;
}

std::vector<std::pair<int, double> > Graph::get_neighbor(int idx) {
    return adjList[idx];
}

std::vector<std::pair<Node, kostMinEdge> > Graph::get_neighbor_KostMin(int idx) {
    return adjListKostMin[idx];
}

std::vector<std::vector<std::pair<int, double> > > Graph::get_adjList() {
    return adjList;
}

std::vector<std::vector<std::pair<Node, kostMinEdge> > > Graph::get_adjListKostMin() {
    return adjListKostMin;
}

std::vector<Edge> Graph::get_EdgeList() {
    return edgeList;
}

std::vector<Node> Graph::get_NodeList() {
    return nodeList;
}

void Graph::adjust_weight(int node1, int node2, double weight) {
    if(!adjMatrix.empty()) {
        adjMatrix[node1][node2] += weight;
    } else if (!adjList.empty()) {
        for (auto& i : adjList[node1]) {
            if (i.first == node2) {
                i.second += weight;
                break;
            }
        }
    } else {
        std::cout << "No AdjList or AdjMatrix initialized!" << std::endl;
    }
}


void Graph::adjust_Kapa(int node1, int node2, double weight) {
    for (auto& i : adjListKostMin[node1]) {
        if (i.second.node2.node1 == node2) {
            i.second.kapazitaet += weight;
            break;
        }
    }
}

void Graph::adjust_Fluss(int node1, int node2, double weight) {
    for (auto& i : adjListKostMin[node1]) {
        if (i.second.node2.node1 == node2) {
            i.second.fluss+= weight;
            break;
        }
    }
}

double Graph::get_Weight(int node1, int node2) {
    if(!adjMatrix.empty()) {
        return adjMatrix[node1][node2];
    } else if (!adjList.empty()) {
        for(auto i : adjList[node1]) {
            if(i.first == node2) {
                return i.second;
            }
        }
    } else {
        std::cout << "No AdjList or AdjMatrix initialized!" << std::endl;
        return 0;
    }
}

// Get the capacity of the edge between node1 and node2
double Graph::get_Kapa(int node1, int node2) {
    for(auto& i : adjListKostMin[node1]) {
        if(i.second.node2.node1 == node2) {
            return i.second.kapazitaet;
        }
    }
}

double Graph::get_Kosten(int node1, int node2) {
    for(auto& i : adjListKostMin[node1]) {
        if(i.second.node2.node1 == node2) {
            return i.second.kosten;
        }
    }
}

double Graph::get_Fluss(int node1, int node2) {
    for(auto& i : adjListKostMin[node1]) {
        if(i.second.node2.node1 == node2) {
            return i.second.fluss;
        }
    }
}

void Graph::addResidualEdges() {
    for (int i = 0; i < V; i++) {
        for (auto& neighbour : get_neighbor(i)) {
            addEdgeToList(neighbour.first, i, 0.0);
        }
    }    
}

void Graph::addResidualEdgesMaxFlow(){
    std::vector<int> x(adjListKostMin.size());
    for (int i = 0; i < adjListKostMin.size(); i++) {
        x[i] = adjListKostMin[i].size();
    }
    for (int i = 0; i < adjListKostMin.size(); i++) {
        for (int j = 0; j < x[i]; j++) {
            Node node1 = {adjListKostMin[i][j].first.node1, adjListKostMin[i][j].first.balance};
            Node node2 = {adjListKostMin[i][j].second.node2.node1, adjListKostMin[i][j].second.node2.balance};
            kostMinEdge edge = {node1, -adjListKostMin[i][j].second.kosten, 0.0, 0.0};
            addEdgeToListKostMin(node2.node1, node2, edge);
        }
    }
}

// Add residual edges to the graph
void Graph::addResidualEdgesCostKapa() {
    std::vector<int> x(adjListKostMin.size());
    for (int i = 0; i < adjListKostMin.size(); i++) {
        x[i] = adjListKostMin[i].size();
    }
    for (int i = 0; i < adjListKostMin.size(); i++) {
        for (int j = 0; j < x[i]; j++) {
            Node node1 = {adjListKostMin[i][j].first.node1, adjListKostMin[i][j].first.balance};
            Node node2 = {adjListKostMin[i][j].second.node2.node1, adjListKostMin[i][j].second.node2.balance};
            if (adjListKostMin[i][j].second.fluss > 0) {
                kostMinEdge edge = {node1, -adjListKostMin[i][j].second.kosten, 0.0, adjListKostMin[i][j].second.fluss};
                addEdgeToListKostMin(node2.node1, node2, edge);
            }
            if ((adjListKostMin[i][j].second.kapazitaet - adjListKostMin[i][j].second.fluss) == 0) {
                adjListKostMin[i].erase(adjListKostMin[i].begin()+j);
                x[i]--;
                j--;
            } else {
                adjust_Kapa(node1.node1, node2.node1, -adjListKostMin[i][j].second.fluss);
            }
        }
    }
}

void Graph::set_NodeList(int node, int balance) {
    nodeList[node].balance = balance;
}

void Graph::adjListKostMin_resize(int size) {
    adjListKostMin.resize(size);
}

void Graph::removeNodeKostMin() {
    adjListKostMin.pop_back();
    adjListKostMin.pop_back();
    int node1 = V;
    int node2 = V+1;
    for (int i = 0; i < adjListKostMin.size(); i++) {
        for (int j = 0; j < adjListKostMin[i].size(); j++) {
            if (adjListKostMin[i][j].second.node2.node1 == node1) {
                adjListKostMin[i].erase(adjListKostMin[i].begin() + j);
            }
            if (adjListKostMin[i][j].second.node2.node1 == node2) {
                adjListKostMin[i].erase(adjListKostMin[i].begin() + j);
            }
        }
    }
}

double Graph::get_GesamtKosten() {
    double gesamtKosten = 0.0;
    for (int i = 0; i < adjListKostMin.size(); i++) {
        for (int j = 0; j < adjListKostMin[i].size(); j++) {
            gesamtKosten += adjListKostMin[i][j].second.kosten * adjListKostMin[i][j].second.fluss;
        }
    }
    return gesamtKosten;
}

bool Graph::checkEdge(int node1, int node2) {
    for (auto& i : adjListKostMin[node1]) {
        if (i.second.node2.node1 == node2) {
            return true;
        }
    }
    return false;
}