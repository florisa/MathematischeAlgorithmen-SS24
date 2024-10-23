#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include "edge.h"

struct Node {
    int node1;
    double balance;
};

struct kostMinEdge {
    Node node2;
    double kosten;
    double fluss;
    double kapazitaet;
};

class Graph {
public:
    Graph();
    Graph(std::string filename, bool weighted, bool kostMinFluss);
    Graph(const Graph& other);
    Graph& operator=(const Graph& other);
    void addEdge(int u, int v, double weight);
    void setAdjList(bool gerichtet);
    void setAdjMatrix(bool gerichtet);
    void addEdgeToList(int u, int v, double weight);
    void addEdgeToListKostMin(int u, Node node, kostMinEdge edge);
    void printAdjList();
    void printAdjListKostMin();
    void printAdjMatrix();
    int get_V();
    std::vector<std::pair<int, double> > get_neighbor(int idx);
    std::vector<std::pair<Node, kostMinEdge> > get_neighbor_KostMin(int idx);
    std::vector<std::vector<std::pair<int, double> > > get_adjList();
    std::vector<std::vector<std::pair<Node, kostMinEdge> > > get_adjListKostMin();
    std::vector<Edge> get_EdgeList();
    std::vector<Node> get_NodeList();
    double get_Weight(int node1, int node2);
    double get_Kapa(int node1, int node2);
    void adjust_weight(int node1, int node2, double weight);
    void adjust_Kapa(int node1, int node2, double weight);
    void adjust_Fluss(int node1, int node2, double weight);
    void addResidualEdges();
    void addResidualEdgesMaxFlow();
    void addResidualEdgesCostKapa();
    void adjListKostMin_resize(int size);
    void removeNodeKostMin();
    void set_NodeList(int node, int balance);
    bool checkEdge(int node1, int node2);
    double get_GesamtKosten();
    double get_Kosten(int node1, int node2);
    double get_Fluss(int node1, int node2);
private:
    int V;
    std::vector<std::vector<std::pair<Node, kostMinEdge> > > adjListKostMin;
    std::vector<std::vector<std::pair<int, double> > > adjList;
    std::vector<std::vector<double> > adjMatrix;
    std::vector<Edge> edgeList;
    std::vector<Node> nodeList;
};

#endif // GRAPH_H