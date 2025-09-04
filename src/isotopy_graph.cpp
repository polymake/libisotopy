#include "isotopy_graph.h"
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <unordered_map>

namespace Isotopy {

Graph::Graph(int delta, const std::vector<bool>&sign_vector, const std::set<std::set<int>>& triangles) {
  std::set<std::pair<int, int>> edges;
  for (const auto& triangle : triangles) {
    if (triangle.size() != 3) {
      throw std::invalid_argument("Each triangle must have exactly 3 vertices.");
    }
    std::vector<int> verts(triangle.begin(), triangle.end());
    edges.insert({verts[0], verts[1]});
    edges.insert({verts[1], verts[2]});
    edges.insert({verts[2], verts[0]});
  }
  *this = Graph(delta, sign_vector, edges);
}


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

  // Create copies to avoid modifying the original data
  auto nbs_copy = nbs;
  auto sides_copy = sides;
  auto component_copy = component;
  
  
  size_t max_iterations = 5; // Prevent infinite loops

  for (size_t iter = 0; iter < max_iterations; ++iter) {
  for (size_t i = 0; i < side_points.size(); ++i) {
    int antipode_i = (i + 2 * delta) % N;
    //First check if the antipode is in the same component
    //If this is the case return the component of the current side point
    if (nbs_copy[i].find(antipode_i) != nbs_copy[i].end()) {
      root = component_copy[side_points[i]];
      return root;
    }

    //Check if the sidepoint is connected to antipodel sides, aka sides 0 and 2 or sides 1 and 3
    if ((sides_copy[i][0] && sides_copy[i][2]) || (sides_copy[i][1] && sides_copy[i][3]))  continue; // Do nothing

    // Check that i together with nbs[i] forms a connected sequence of integers (e.g., {3,4,5,6}),
    // considering indices as circular (0 is next to 4*delta-1)
    std::set<int> all_indices = nbs_copy[i];
    all_indices.insert(i);
    int min_idx = *std::min_element(all_indices.begin(), all_indices.end());
    int max_idx = *std::max_element(all_indices.begin(), all_indices.end());
    bool is_connected = false;
    if ((max_idx - min_idx + 1) == (int)all_indices.size()) {
      is_connected = true; // linear connected
    } else {
      //THIS IS AN ISSUE
      //PLEASE FIX
      //THis ignores wrap around cases, e,g {29,30,31,0,1} for delta = 8
      //NOTE: I believe this is not an issue since it will just merge the opposite side points.
      is_connected = false;
    }
    /*
    if (i < 5) {
      std::cout << "Iteration: " << iter << "\n";
      std::cout << "i: " << i << " nbs: ";
      for (int nb : nbs_copy[i]) {
        std::cout << nb << " ";
      }
      std::cout << " is_connected: " << is_connected << "\n";
    }
    */

    if (!is_connected) continue; // Do nothing 

    int j = (i-1) % N;
    //merge the nbs of i and j
    std::set<int> merged_nbs;
    merged_nbs.insert(nbs_copy[i].begin(), nbs_copy[i].end());
    merged_nbs.insert(nbs_copy[j].begin(), nbs_copy[j].end());
    merged_nbs.insert(i);
    merged_nbs.insert(j);

    //Update sides of j to include sides of i
    std::vector<bool> new_sides(4, false);

    for (size_t k = 0; k < 4; ++k) {
      new_sides[k] = sides_copy[i][k] || sides_copy[j][k];
    }
    //Change the component of side point i to that of j
    component_copy[side_points[i]] = component_copy[side_points[j]];
    //Check if antipode of j is in nbs[j] or the antipode of i is in nbs[j]
    int antipode_j = (j + 2 * delta) % N;
    if (nbs_copy[j].find(antipode_j) != nbs_copy[j].end() || nbs_copy[j].find(antipode_i) != nbs_copy[j].end()) {
      root = component_copy[side_points[j]];
      return root;
    }
    //Update all elements in merged_nbs to have the same nbs and sides
    for (int idx : merged_nbs) {
      nbs_copy[idx] = merged_nbs;
      sides_copy[idx] = new_sides;
      component_copy[side_points[idx]] = component_copy[side_points[j]];
    }

    //Merge antipode_i with antipode_j
    std::set<int> merged_nbs_anti;
    merged_nbs_anti = nbs_copy[antipode_j];
    merged_nbs_anti.insert(antipode_i);
     //Update sides of antipode_j to include sides of antipode_i
    std::vector<bool> new_sides_anti(4, false);

    for (size_t k = 0; k < 4; ++k) {
      new_sides_anti[k] = sides[antipode_i][k] || sides_copy[antipode_j][k];
    }

    component_copy[side_points[antipode_i]] = component_copy[side_points[antipode_j]];
    
    for (int idx : merged_nbs_anti) {
      nbs_copy[idx] = merged_nbs_anti;
      sides_copy[idx] = new_sides_anti;
    }

  }
  }
  std::cerr << "Warning: No isotopy root found after maximum iterations.\n";
  return -1; // No isotopy root found
}

void Graph::calculate_regions() {
  if (root == -1) {
    root = isotopy_root();
  }

  // Build component_antipod_adjacency
  std::vector<std::set<int>> component_antipode_adjacency(component_adjacency.size());
  for (size_t i = 0; i < side_points.size(); ++i) {
    int comp = component[side_points[i]];
    int antipodal_comp = component[side_points[(i + 2 * delta) % (4 * delta)]];
    component_antipode_adjacency[comp].insert(antipodal_comp);
  }

  region.resize(component_adjacency.size(), -1);
  int curr = 0;
  for (size_t i = 0; i < region.size(); ++i) {
    if (region[i] == -1) {
      std::set<int> stack;
      stack.insert(i);
      while (!stack.empty()) {
        int current = *stack.begin();
        stack.erase(stack.begin());
        region[current] = curr;
        for (int neighbor : component_antipode_adjacency[current]) {
          if (region[neighbor] == -1) {
            stack.insert(neighbor);
          }
        }
      }
      curr++;
    }
  }


  region_adjacency.resize(curr, std::set<int>());
  for (size_t i = 0; i < component_adjacency.size(); ++i) {
    for (int adj_comp : component_adjacency[i]) {
      region_adjacency[region[i]].insert(region[adj_comp]);
    }
  }
  root_region = region[root];
}

void Graph::isotopy_type() {
  if (component.empty()) {
    connected_components();
  }
  if (root == -1) {
    pre_isotopy_root();
    root = isotopy_root();
    if (root == -1) {
      throw std::runtime_error("No isotopy root found");
    }
  }
  if (region_adjacency.empty()) {
    calculate_regions();
  }



  region_sign.resize(region_adjacency.size(), false);

  std::set<int> stack;
  stack.insert(root_region);
  region_sign[root_region] = false; 
  std::set<int> visited;

  while (!stack.empty()) {
    int current = *stack.begin();
    stack.erase(stack.begin());
    visited.insert(current);
    if (region_sign[current]) {
      p_regions++;
    } else {
      n_regions++;
    }

    for (int neighbor : region_adjacency[current]) {
      if (visited.find(neighbor) == visited.end()) {
        region_sign[neighbor] = !region_sign[current];
        edges.insert({current, neighbor});
        stack.insert(neighbor);
      }
    }
  }
}

}
