#include "isotopy_graph.h"
#include <cassert>
#include <stdexcept>
#include <algorithm>

namespace Isotopy {

Graph::Graph(int delta, const std::vector<bool>& sign, const std::set<std::pair<int, int>>& edges)
    : delta(delta) {
    std::size_t nverts = (delta + 1) * (delta + 2) / 2;
    assert(sign.size() == nverts && "sign vector length does not match number of vertices");

    for (const auto& edge : edges) {
        assert(edge.first < static_cast<int>(nverts) && edge.second < static_cast<int>(nverts) && "edge vertex index out of range");
        if (sign[edge.first] != sign[edge.second]) continue;
        adjacency[edge.first].insert(edge.second);
        adjacency[edge.second].insert(edge.first);
    }

}


}



