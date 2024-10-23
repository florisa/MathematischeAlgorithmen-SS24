#include "graph.h"
#include "edge.h"
#include "algo.h"
#include <chrono>
 
int main() {
    
    Graph g1("Kostenminimal1.txt", true, true);
    Graph g2("Kostenminimal2.txt", true, true);
    Graph g3("Kostenminimal3.txt", true, true);
    Graph g4("Kostenminimal4.txt", true, true);
    Graph g5("Kostenminimal_gross1.txt", true, true);
    Graph g6("Kostenminimal_gross2.txt", true, true);
    Graph g7("Kostenminimal_gross3.txt", true, true);
    std::cout << "All files read!" << std::endl;

    std::vector<Graph> graphs;
    graphs.push_back(g1);
    graphs.push_back(g2);
    graphs.push_back(g3);
    graphs.push_back(g4);
    graphs.push_back(g5);
    graphs.push_back(g6);
    graphs.push_back(g7);


    for(Graph& graph : graphs) {
        //double minFluss= successive_shortest_Path(graph);
        double minFluss2 = cycle_canceling(graph);
        //std::cout << "Kostenminimaler Fluss ist " << minFluss << std::endl << std::endl << std::endl;
        std::cout << "Kostenminimaler Fluss2 ist " << minFluss2 << std::endl << std::endl << std::endl;
    }

    return 0;
}