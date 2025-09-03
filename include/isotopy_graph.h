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
    std::vector<std::vector<bool>> sides; // sides[i][j] is true if sidepoint i is on side j (0: top-right, 1: bottom-right, 2: bottom-left, 3: top-left)
    std::vector<std::set<int>> nbs; // nbs[i] gives the set of other sidepoints connected to sidepoint i via its component

    // Information about connected components
    int root = -1; // root component index
    std::vector<int> component; // component[i] gives the component index of vertex i
    std::vector<std::set<int>> component_adjacency; // component_adjacency[c] gives the set of components adjacent to component c

    
    Graph() = default;
    Graph(int delta, const std::vector<bool>& sign_vector, const std::set<std::pair<int, int>>& edges);

    void connected_components();
    void pre_isotopy_root();
    int isotopy_root();
    
};



}
