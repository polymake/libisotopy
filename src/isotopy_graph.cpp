#include "isotopy_graph.h"
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <unordered_map>

namespace Isotopy {

Graph::Graph(int delta, const std::vector<bool>&sign_vector, const std::set<std::pair<int, int>>& edges)
: delta(delta) {
  std::size_t nverts = (delta + 1) * (delta + 2) / 2;
  assert(sign_vector.size() == nverts && "sign vector length does not match number of vertices");

  size_t nverts_triang = (delta + 1) * (delta + 2) / 2;
  size_t nvert_formula = 4 * nverts_triang - 4 * (delta + 1) +1;

  std::map<int,int> vertex_map;

  sign.resize(nvert_formula);
  side_points.resize(4 * delta);
  sides.resize(4 * delta, std::vector<bool>(4, false)); 

  int side_point_count = 0;
  int global_count = 0;
  int count = 0;
  //First Quadrant
  for (int limit = delta; limit >= 0; --limit) {
    for (int i = 0; i <= limit; ++i) {
      if (i == limit) {
        side_points[side_point_count] = count;
        if (limit == 0) {
          sides[side_point_count][0] = true;
          sides[side_point_count][1] = true;
        } else if (limit == delta) {
          sides[side_point_count][0] = true;
          sides[side_point_count][3] = true;
        } else {
          sides[side_point_count][0] = true;
        }
        side_point_count++;

      }

      vertex_map[count] = count;
      sign[count] =sign_vector[count];
      count++;
      global_count++;
    }
  }
  //Second Quadrant
  count = 0;
  for (int limit = delta; limit >= 0; --limit) {
    for (int i = 0; i <= limit; ++i) {
      if (i == 0) {
        vertex_map[count+nverts] = count;
      } else {
        if (i == limit) {
          side_points[side_point_count +  limit -1 ] = global_count;
          if (limit == delta) {
            sides[side_point_count + limit -1][1] = true;
            sides[side_point_count + limit -1][2] = true;
          } else {
            sides[side_point_count + limit -1][1] = true;
          }
        }
        vertex_map[count+nverts] = global_count;
        if (i % 2 == 0) {
          sign[global_count] = sign_vector[count];
        } else {
          sign[global_count] = !sign_vector[count];
        }
        global_count++;
      }
      count++;
    }
  }
  side_point_count += delta;
  
  //Third and Fourth Quadrants
  count = 0;
  for (int limit = delta; limit >= 0; --limit) {
    for (int i = 0; i <= limit; ++i) {
      if (limit == delta) {
        vertex_map[count+2*nverts] = vertex_map[count+nverts];
        vertex_map[count+3*nverts] = vertex_map[count];
      } else if (i == 0) {
        vertex_map[count+2*nverts] = global_count;
        vertex_map[count+3*nverts] = global_count;
        if (limit % 2 == 0) {
          sign[global_count] = sign_vector[count];
        } else {
          sign[global_count] = !sign_vector[count];
        }
        if (limit == 0) {
          side_points[side_point_count] = global_count;
          sides[side_point_count][2] = true;
          sides[side_point_count][3] = true;
          side_point_count++;
        }
        global_count++;
      } else {
        vertex_map[count+2*nverts] = global_count;
        if (i % 2 == 0) {
          sign[global_count] = sign_vector[count];
        } else {
          sign[global_count] = !sign_vector[count];
        }
        if (limit % 2 == 0) {
          sign[global_count] = sign[global_count];
        } else {
          sign[global_count] = !sign[global_count];
        }
        //Save new side points 
        if (i == limit) {
          side_points[side_point_count] = global_count;
          sides[side_point_count][2] = true;
          side_points[side_point_count + 2*limit] = global_count+1;
          sides[side_point_count + 2*limit][3] = true;
          side_point_count++;
        }



        global_count++;
        vertex_map[count+3*nverts] = global_count;
        if (limit % 2 == 0) {
          sign[global_count] = sign_vector[count];
        } else {
          sign[global_count] = !sign_vector[count];
        }
        global_count++;
      }
      count++;
    }
  }
  /*
  std::cout << "Vertex Map:\n";
  for (const auto& [key, value] : vertex_map) {
    std::cout << "Original: " << key << " -> New: " << value << "\n";
  }
  */

  adjacency.resize(global_count);
  for (const auto& edge : edges) {
    assert(edge.first < static_cast<int>(nverts) && edge.second < static_cast<int>(nverts) && "edge vertex index out of range");
    for (size_t k = 0; k < 4; ++k) {
      int v1 = vertex_map[edge.first + k * nverts];
      int v2 = vertex_map[edge.second + k * nverts];

      adjacency[v1].insert(v2);
      adjacency[v2].insert(v1);
    }
  }
}

void Graph::connected_components() {
  std::set<int> visited;
  int component_count = 0;

  component.resize(adjacency.size(), -1); // Ensure enough space for all vertex indices

  for (size_t vertex = 0; vertex < adjacency.size(); ++vertex) {
    if (visited.find(vertex) == visited.end()) {
      std::set<int> stack;
      stack.insert(vertex);
      std::set<int> temp_adjacency =  std::set<int>();

      while (!stack.empty()) {
        int current = *stack.begin();
        stack.erase(stack.begin());
        visited.insert(current);
        component[current] = component_count;

        for (int neighbor : adjacency[current]) {
          if (visited.find(neighbor) == visited.end()) {
            if (sign[neighbor] != sign[current]) continue; // Only traverse edges with same sign
              stack.insert(neighbor);
          } else {
            if (sign[neighbor] != sign[current]){
              auto c2 = component[neighbor];
              temp_adjacency.insert(c2);
              component_adjacency[c2].insert(component_count);

            }
          }
        }
      }
      component_adjacency.push_back(temp_adjacency);
      component_count++;
    }
  }

}

void Graph::pre_isotopy_root() {
  /*
  We need to compute: 
  nbs: A vector of sets of other sidepoints a point is connected to via its component
  */
  nbs.resize(side_points.size());
  std::unordered_map<int, std::vector<int>> comp_to_indices;
  for (size_t i = 0; i < side_points.size(); ++i) {
    int vertex = side_points[i];
    int comp = component[vertex];
    comp_to_indices[comp].push_back(i);
  }

  for (const auto& pair : comp_to_indices) {
    const std::vector<int>& indices = pair.second;
    for (int idx : indices) {
      for (int other_idx : indices) {
        if (idx != other_idx) {
          nbs[idx].insert(other_idx);
        }
      }
    }
  }
  //int antipode = side_points[(i + 2 * delta) % (4 * delta)];
}


int Graph::isotopy_root() {

  int N = 4 * delta; // total number of indices

  for (size_t i = 0; i < side_points.size(); ++i) {
    int antipode_i = (i + 2 * delta) % N;
    //First check if the antipode is in the same component
    //If this is the case return the component of the current side point
    if (nbs[i].find(antipode_i) != nbs[i].end()) {
      root = component[side_points[i]];
      return root;
    }

    //Check if the sidepoint is connected to antipodel sides, aka sides 0 and 2 or sides 1 and 3
    if ((sides[i][0] && sides[i][2]) || (sides[i][1] && sides[i][3]))  continue; // Do nothing

    // Check that i together with nbs[i] forms a connected sequence of integers (e.g., {3,4,5,6}),
    // considering indices as circular (0 is next to 4*delta-1)
    std::set<int> all_indices = nbs[i];
    all_indices.insert(i);
    int min_idx = *std::min_element(all_indices.begin(), all_indices.end());
    int max_idx = *std::max_element(all_indices.begin(), all_indices.end());
    bool is_connected = false;
    if ((max_idx - min_idx + 1) == (int)all_indices.size()) {
      is_connected = true; // linear connected
    } else if ((int)all_indices.size() == N - (max_idx - min_idx - 1)) {
      // circular connected: indices wrap around
      // e.g., {N-2, N-1, 0, 1}
      is_connected = true;
    }

    if (!is_connected) continue; // Do nothing 

    int j = (i-1) % N;
    //merge the nbs of i and j
    nbs[j].insert(nbs[i].begin(), nbs[i].end());
    nbs[i] = nbs[j];
    //Update sides of j to include sides of i
    for (size_t k = 0; k < 4; ++k) {
      sides[j][k] = sides[j][k] || sides[i][k];
      sides[i][k] = sides[j][k];
    }
    //Change the component of side point i to that of j
    component[side_points[i]] = component[side_points[j]];
    //Check if antipode of j is in nbs[j] or the antipode of i is in nbs[j]
    int antipode_j = (j + 2 * delta) % N;
    if (nbs[j].find(antipode_j) != nbs[j].end() || nbs[j].find(antipode_i) != nbs[j].end()) {
      root = component[side_points[j]];
      return root;
    }

    //Merge antipode_i with antipode_j
    nbs[antipode_j].insert(nbs[antipode_i].begin(), nbs[antipode_i].end());
    nbs[antipode_i] = nbs[antipode_j];
    //Update sides of antipode_j to include sides of antipode_i
    for (size_t k = 0; k < 4; ++k) {
      sides[antipode_j][k] = sides[antipode_j][k] || sides[antipode_i][k];
      sides[antipode_i][k] = sides[antipode_j][k];
    }
  }
  return -1; // No isotopy root found
}

}
