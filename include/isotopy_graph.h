#pragma once

#include <vector>
#include <set>
#include <map>
//#include <algorithm>
//#include <cstddef>

namespace Isotopy {

struct Graph {
    int delta;
    std::map<int, std::set<int>> adjacency;

    
    Graph() = default;
    Graph(int delta, const std::vector<bool>& sign, const std::set<std::pair<int, int>>& edges);
    
};

}
