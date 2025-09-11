#include "isotopy_graph.h"
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <functional>

//For the Utils
#include <cstdint>

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
  if (component.empty()) {
    connected_components();
  }


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
  if (nbs.empty()) {
    pre_isotopy_root();
  }

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

    int j = ((i+N)-1) % N;
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
}

void Graph::isotopy_type() {
  if (root == -1) {
    isotopy_root();
  }

  if (region_adjacency.empty()) {
    calculate_regions();
  }

  root_region = region[root];

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

std::string Graph::viro_notation(bool unicode) {
  if (region_adjacency.empty()) {
    isotopy_type();
  }
  return Isotopy::viro_notation(region[root], region_adjacency, unicode);
}

std::string viro_notation(int root_region, const std::vector<std::set<int>>& region_adjacency, bool unicode) {
  const std::string open_delim = unicode ? "\u27E8" : "<";
  const std::string close_delim = unicode ? "\u27E9" : ">";
  const std::string sep = unicode ? "\u2294" : "v";

  std::function<std::string(int, std::set<int>&)> dfs =
    [&](int curr_region, std::set<int>& visited) -> std::string {
      visited.insert(curr_region);
      int leaf_count = 0;
      std::vector<int> non_leaf_children;
      for (const auto& neighbor : region_adjacency[curr_region]) {
        if (!visited.count(neighbor)) {
          bool is_leaf = true;
          for (const auto& nn : region_adjacency[neighbor]) {
            if (!visited.count(nn) && nn != curr_region) {
              is_leaf = false;
              break;
            }
          }
          if (is_leaf) {
            leaf_count++;
            visited.insert(neighbor);
          } else {
            non_leaf_children.push_back(neighbor);
          }
        }
      }
      if (non_leaf_children.empty()) {
        return open_delim + std::to_string(leaf_count) + close_delim;
      } else if (leaf_count == 0 && non_leaf_children.size() == 1) {
        return open_delim + "1" + dfs(non_leaf_children[0], visited) + close_delim;
      } else {
        std::vector<std::string> child_types;
        for (const auto& child : non_leaf_children) {
          child_types.push_back(dfs(child, visited));
        }
        std::map<std::string, int> type_counts;
        for (const auto& t : child_types) {
          type_counts[t]++;
        }
        std::vector<std::pair<int, std::string>> count_type_pairs;
        for (auto it = type_counts.begin(); it != type_counts.end(); ++it) {
          count_type_pairs.emplace_back(it->second, it->first);
        }
        std::sort(count_type_pairs.begin(), count_type_pairs.end());
        std::vector<std::string> grouped_types;
        for (const auto& pair : count_type_pairs) {
          grouped_types.push_back(std::to_string(pair.first) + pair.second);
        }
        std::string result = open_delim;
        if (leaf_count > 0) {
          result += std::to_string(leaf_count) + sep;
        }
        for (size_t i = 0; i < grouped_types.size(); ++i) {
          if (i > 0) result += sep;
          result += grouped_types[i];
        }
        result += close_delim;
        return result;
      }
    };

  std::set<int> visited;
  return dfs(root_region, visited);
}
}

namespace Utils {

//The website enumerates the points from [0,delta]->0, [0,delta-1]->1, ..., [0,1]->delta-1, [0,0]->delta
//whereas this does it the other way around

std::vector<std::pair<int,int>> lib2pt(int delta) {
  std::vector<std::pair<int,int>> vertex_map((delta + 1) * (delta + 2) / 2);
  size_t count = 0;
  for (int limit = delta; limit >= 0; --limit) {
    for (int i = 0; i <= limit; ++i) {
      vertex_map[count] = {delta - limit, i};
      count++;
    }
  }
  return vertex_map;
}
std::vector<std::pair<int,int>> web2pt(int delta) {
std::vector<std::pair<int,int>> vertex_map((delta + 1) * (delta + 2) / 2);
size_t count = 0;
for (int x = 0; x <= delta; ++x) {
  for (int y = delta - x; y <= delta; ++y) {
    vertex_map[count] = {x,y};
    count++;
  }
}
return vertex_map;
}
std::map<std::pair<int,int>,int> pt2web(int delta) {
  // The order in web2pt is: for x in [0, delta], for y in [delta-x, delta]
  // So for a given (x, y), the index is:
  // sum_{x'=0}^{x-1} (x'+1) + (y - (delta - x))
  
  std::map<std::pair<int,int>,int> pt2web;
  int count = 0;
  for (int limit = delta; limit >= 0; --limit) {
    for (int y = limit; y>=0; --y) {
      std::cout << y << "\n";
      int x = delta - limit;
      pt2web[std::make_pair(x,y)] = count;
      count++;
    }
  }
  return pt2web;
}

std::vector<int> lib2web(int delta) {
  auto lib_pts = lib2pt(delta);
  std::cout << "hello\n";
  auto pt_2web = pt2web(delta);
  std::cout <<"wat\n";
  std::vector<int> mapping(lib_pts.size());
  for (size_t i = 0; i < lib_pts.size(); ++i) {
    std::cout << "Creating lin_pts\n";
    int x = lib_pts[i].first;
    int y = lib_pts[i].second;
    mapping[i] = pt_2web[std::make_pair(x,y)];
  }
  return mapping;
}

std::vector<std::vector<int>>  adjacency_matrix(int delta, const std::set<std::pair<int, int>>& edges, bool web_format) {
  size_t n = (delta + 1) * (delta + 2) / 2;
  std::vector<std::vector<int>> matrix(n, std::vector<int>(n, 0));

  for (const auto& edge : edges) {
    int u = edge.first;
    int v = edge.second;
    if (u < 0 || u >= static_cast<int>(n) || v < 0 || v >= static_cast<int>(n)) {
      throw std::out_of_range("Edge vertex index out of range");
    }
    matrix[u][v] = 1;
    matrix[v][u] = 1; // Undirected graph
  }

  if (web_format) {
    std::vector<int> mapping = lib2web(delta);
    std::cout << mapping.size() << "\n";
    std::cout << n << "\n";
    std::vector<std::vector<int>> web_matrix(n, std::vector<int>(n, 0));
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        std::cout << "hello " << i << " " << j << "\n";
        std::cout << mapping[i] << " , " << mapping[j] << "\n";
        web_matrix[mapping[i]][mapping[j]] = matrix[i][j];
      }
    }
    return web_matrix;
  }

  return matrix;
}

std::string to_graph6(const std::vector<std::vector<int>>& adjacency_matrix) {
  int n = adjacency_matrix.size();
  // Use std::vector<uint8_t> for better cache locality and bitwise ops
  std::vector<uint8_t> edge_bits;
  edge_bits.reserve(n * (n - 1) / 2);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      edge_bits.push_back(adjacency_matrix[i][j] == 1 ? 1 : 0);
    }
  }
  // Pad right so the length is a multiple of 6
  int padding = (6 - edge_bits.size() % 6) % 6;
  edge_bits.insert(edge_bits.end(), padding, 0);

  std::string graph6;
  for (size_t i = 0; i < edge_bits.size(); i += 6) {
    uint8_t value = 0;
    for (int j = 0; j < 6; ++j) {
      value = (value << 1) | edge_bits[i + j];
    }
    graph6 += static_cast<char>(value + 63);
  }
  return graph6;
}

  
std::string shorthand(int delta, std::vector<bool> sign_vector, std::set<std::pair<int, int>> edges) {
  // int nedges =  3*delta+3/2 * (delta*delta-delta); // Unused
  // int nfaces = delta*delta; // Unused

  int nverts = (delta + 1) * (delta + 2) / 2; // Number of vertices in the triangulation
  int length_of_adjacency = (nverts * (nverts - 1)) / 2; // Number of edges in a complete graph
  int alt_padding = (6 - length_of_adjacency % 6) % 6; // Pad to multiple of 6
  int length_of_graph6 = (length_of_adjacency + alt_padding) / 6;
  int padding_sign_str = (6 - nverts % 6) % 6; // Pad to multiple of 6
  int length_of_sign_str = (nverts + padding_sign_str) / 6;

  auto adj_matrix = adjacency_matrix(delta, edges);

  //If web format convert the sign vector
  /*
  if (web_format) {
    std::vector<bool> sign_vector_web(nverts);
    std::vector<int> mapping = lib2web(delta);
    for (int i = 0; i < nverts; ++i) {
      sign_vector_web[mapping[i]] = sign_vector[i];
    }
    sign_vector = sign_vector_web;
  }
  */

  //Convert the sign string.
  size_t padding = (6 - sign_vector.size() % 6) % 6;
  std::vector<bool> sign_vector_padded(sign_vector.size() + padding, 0);
  if (static_cast<size_t>(length_of_sign_str) != (sign_vector_padded.size() / 6)) {
    std::cerr << "Error: Length of sign string does not match expected length." << std::endl;
    return std::string();
  }
  for (size_t i = 0; i < sign_vector.size(); ++i) {
    sign_vector_padded[i] = sign_vector[i];
  }
  std::string sign_short;
  for (size_t i = 0; i < sign_vector_padded.size(); i += 6) {
    uint8_t value = 0;
    for (int j = 0; j < 6; ++j) {
      value = (value << 1) | sign_vector_padded[i + j];
    }
    sign_short += static_cast<char>(value + 63);
  }

  std::string g6 = to_graph6(adj_matrix);

  //As the leading character we have the delta+63 as ascii
  std::string shorthand_string = std::string(1, static_cast<char>(delta + 63)) + g6 + sign_short;
  if (shorthand_string.length() != static_cast<std::string::size_type>(1 + length_of_graph6 + length_of_sign_str)) {
    std::cerr << "Error: Shorthand string length does not match expected length." << std::endl;
  }

  return shorthand_string;

}

std::string shorthand(int delta, std::vector<bool> sign_vector, std::set<std::set<int>> triangles) {
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
  return shorthand(delta, sign_vector, edges);
}

std::string web_shorthand(int delta, std::vector<bool> sign_vector, std::set<std::pair<int, int>> edges) {
  // int nedges =  3*delta+3/2 * (delta*delta-delta); // Unused
  // int nfaces = delta*delta; // Unused

  int nverts = (delta + 1) * (delta + 2) / 2; // Number of vertices in the triangulation
  int length_of_adjacency = (nverts * (nverts - 1)) / 2; // Number of edges in a complete graph
  int alt_padding = (6 - length_of_adjacency % 6) % 6; // Pad to multiple of 6
  int length_of_graph6 = (length_of_adjacency + alt_padding) / 6;
  int padding_sign_str = (6 - nverts % 6) % 6; // Pad to multiple of 6
  int length_of_sign_str = (nverts + padding_sign_str) / 6;

  auto adj_matrix = adjacency_matrix(delta, edges,true);
  std::cout<< "Survived adjaceceny";
  //If web format convert the sign vector
  std::vector<bool> sign_vector_web(nverts);
  std::vector<int> mapping = lib2web(delta);
  for (int i = 0; i < nverts; ++i) {
    sign_vector_web[mapping[i]] = sign_vector[i];
  }
  sign_vector = sign_vector_web;

  //Convert the sign string.
  size_t padding = (6 - sign_vector.size() % 6) % 6;
  std::vector<bool> sign_vector_padded(sign_vector.size() + padding, 0);
  if (static_cast<size_t>(length_of_sign_str) != (sign_vector_padded.size() / 6)) {
    std::cerr << "Error: Length of sign string does not match expected length." << std::endl;
    return std::string();
  }
  for (size_t i = 0; i < sign_vector.size(); ++i) {
    sign_vector_padded[i] = sign_vector[i];
  }
  std::string sign_short;
  for (size_t i = 0; i < sign_vector_padded.size(); i += 6) {
    uint8_t value = 0;
    for (int j = 0; j < 6; ++j) {
      value = (value << 1) | sign_vector_padded[i + j];
    }
    sign_short += static_cast<char>(value + 63);
  }

  std::string g6 = to_graph6(adj_matrix);

  //As the leading character we have the delta+63 as ascii
  std::string shorthand_string = std::string(1, static_cast<char>(delta + 63)) + g6 + sign_short;
  if (shorthand_string.length() != static_cast<std::string::size_type>(1 + length_of_graph6 + length_of_sign_str)) {
    std::cerr << "Error: Shorthand string length does not match expected length." << std::endl;
  }

  return shorthand_string;

}

std::string web_shorthand(int delta, std::vector<bool> sign_vector, std::set<std::set<int>> triangles) {
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
  return web_shorthand(delta, sign_vector, edges);
}

/*
std::tuple<int, std::vector<bool>, std::set<std::pair<int, int>>> parse_shorthand(const std::string& shorthand) {

}
*/
}

