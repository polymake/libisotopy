#include "isotopy_graph.h"
#include <cassert>
#include <iostream>

int main() {
    
    int delta = 2;
    std::vector<bool> sign = {true, true, true, true, true, true};
    std::set<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {0, 3}, {3, 4}, {2, 4}, {1,3}, {1,4}, {4,5}, {3,5}};
    std::cout << "length of edges: " << edges.size() << "\n";

    Isotopy::Graph graph = Isotopy::Graph(delta, sign, edges);

    assert(graph.adjacency.size() == 13);
    size_t nverts_triang = (delta + 1) * (delta + 2) / 2;
    size_t nvert_formula = 4 * nverts_triang - 4 * (delta + 1) +1;
    assert(graph.adjacency.size() == nvert_formula);
    
    for (size_t vertex = 0; vertex < graph.adjacency.size(); ++vertex) {
        std::cout << "Vertex " << vertex << ": ";
        for (const auto& neighbor : graph.adjacency[vertex]) {
            std::cout << neighbor << " ";
        }
        std::cout << "\n";
    }

    
    // Adjacency size checks
    assert(graph.adjacency[0].size() == 4);
    assert(graph.adjacency[1].size() == 6);
    assert(graph.adjacency[10].size() == 4);
    assert(graph.adjacency[12].size() == 3);
    assert(graph.adjacency[9].size() == 6);
    
    /*
    std::cout << "Sign vector:\n";
    for (size_t i = 0; i < graph.adjacency.size(); ++i) {
        std::cout << "Vertex " << i << ": Sign " << graph.sign[i] << "\n";
    }
    */
    //Sign vector checks
    assert(graph.sign == std::vector<bool>({true, true, true, true, true, true, false, true, false, false, true, false, true}));


    graph.connected_components();
    
    assert(graph.component.size() == graph.adjacency.size());
    assert(graph.component_adjacency.size() == 3);
    assert(graph.component == std::vector<int>({0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 2, 1, 2}));

    //Print all components
    std::cout << "Number of components: " << graph.component_adjacency.size() << "\n";
    for (size_t c = 0; c < graph.component_adjacency.size(); ++c) {
        std::cout << "Component " << c << ": ";
        for (const auto& neighbor_comp : graph.component_adjacency[c]) {
            std::cout << neighbor_comp << " ";
        }
        std::cout << "\n";
    }
    for (size_t v = 0; v < graph.component.size(); ++v) {
        std::cout << "Vertex " << v << " is in component " << graph.component[v] << "\n";
    }



    return 0;
}
