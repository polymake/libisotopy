#include "isotopy_graph.h"
#include <cassert>
#include <iostream>

int main() {
    
    int delta = 2;
    std::vector<bool> sign = {true, true, true, true, true, true};
    std::set<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {0, 3}, {3, 4}, {2, 4}, {1,3}, {1,4}, {4,5}, {3,5}};
    std::cout << "length of edges: " << edges.size() << "\n";

    Isotopy::Graph graph = Isotopy::Graph(delta, sign, edges);

    graph.adjacency[0].insert(1);
    graph.adjacency[1].insert(0);

    assert(graph.adjacency.size() == 6);
    // 0 should be adjacent to 2 vertices
    assert(graph.adjacency[0].size() == 2);
    assert(graph.adjacency[1].size() == 4);

    std::cout << "Dummy tests passed.\n";
    return 0;
}
