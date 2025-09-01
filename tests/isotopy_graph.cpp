#include "isotopy_graph.h"
#include <cassert>
#include <iostream>

int main() {
    IsotopyGraph graph;
    graph.adjacency.resize(2);
    graph.adjacency[0].insert(1);
    graph.adjacency[1].insert(0);

    assert(graph.adjacency.size() == 2);
    assert(graph.adjacency[0].count(1) == 1);
    assert(graph.adjacency[1].count(0) == 1);

    std::cout << "Dummy tests passed.\n";
    return 0;
}
