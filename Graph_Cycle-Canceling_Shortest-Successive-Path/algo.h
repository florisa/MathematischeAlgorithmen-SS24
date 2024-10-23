#ifndef ALGO_H
#define ALGO_H

#include "graph.h"
#include "edge.h"
#include <stack>
#include <queue>
#include <list>
#include <limits>
#include <algorithm>
#include <iterator>

std::vector<bool> depth_first_search(Graph& graph, int start, std::vector<bool>& visited);

std::vector<int> breadth_first_search(Graph& graph, int start);

std::vector<int> breadth_first_search_Max_Flow(Graph& graph, int start, std::vector<bool>& visited);

std::vector<std::vector<int> > component_bfs(Graph& graph);

int connected_components(Graph& graph);

void prim(Graph& graph, int node);

std::vector<std::vector<int> > prim_tree(Graph& graph, int node);

int find_set(std::vector<int>& parent, int node);

bool union_set(std::vector<int>& parent, std::vector<int>& rank, int node1, int node2);

void kruskal(Graph& graph);

void nearest_neighbour(Graph& graph, int start);

void doppelter_baum(Graph& graph, int start);

void print_matrix(std::vector<std::vector<double> >& matrix);

void alle_hamilton_kreise(Graph& graph);

void generiere_hamilton_kreis(Graph& graph, int current, std::vector<int>& hamilton_kreis, std::vector<bool>& visited, double& hamilton_kreis_weight, std::vector<int>& best_hamilton_kreis, double& best_hamilton_kreis_weight, int& count_hamil_kreise);

void branch_and_bound_hamilton(Graph& graph);

void generiere_hamilton_kreis_b_b(Graph& graph, int current, std::vector<int>& hamilton_kreis, std::vector<bool>& visited, double& hamilton_kreis_weight, std::vector<int>& best_hamilton_kreis, double& best_hamilton_kreis_weight, int& count_hamil_kreise);

void dijkstra(Graph& graph, int startknoten);

void moore_bellman_ford(Graph& graph, int startknoten);

std::vector<int> breadth_first_search_to_target(Graph& graph, int source, int target, double& kleinsterFluss);

std::vector<int> breadth_first_search_to_target_KostMinFlow(Graph& graph, int source, int target, double& kleinsterFluss);

double edmonds_karp(Graph& graph, int source, int target);

double edmonds_karp_KostMinFlowGraph(Graph& graph, int source, int target);

void addSuperSourceAndSink(Graph& graph);

std::vector<int> moore_bellman_ford_negative_cycle(Graph& graph);

std::vector<int> mbf(Graph& graph, int& source, int& sink);

bool graph_Loesbarkeit(Graph& graph, double maxFluss);

void updateFlow(Graph& originalGraph, Graph& graph);

double cycle_canceling(Graph& graph);

bool balancen_ausgelichen(std::vector<Node> nodeList, std::vector<double>& new_balance);

double successive_shortest_Path(Graph& graph);

#endif // ALGO_H