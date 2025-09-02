#pragma once

#include <vector>
#include <set>
#include <map>
//#include <algorithm>
//#include <cstddef>

namespace Isotopy {

struct Graph {
    int delta;
    std::vector<bool> sign;
    std::vector<std::set<int>> adjacency;
    std::vector<int> side_points; 

    // Information about connected components
    std::vector<int> component; // component[i] gives the component index of vertex i
    std::vector<std::set<int>> component_adjacency; // component_adjacency[c] gives the set of components adjacent to component c

    
    Graph() = default;
    Graph(int delta, const std::vector<bool>& sign_vector, const std::set<std::pair<int, int>>& edges);

    void connected_components();
    
};



}
