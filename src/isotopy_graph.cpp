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
  size_t nverts = (delta + 1) * (delta + 2) / 2;
  assert(sign_vector.size() == nverts && "sign vector length does not match number of vertices");

  sign_complete.reserve(4 * nverts);
  edges_complete.reserve(4 * edges.size() + 4*(delta*1));
  for (int i = 0; i < 4; ++i) {
    for (const auto& edge : edges) {
      edges_complete.push_back({edge.first + i * nverts, edge.second + i * nverts});
    }
  }

  //This is wrong
  int offset = 0;
  for (int i = 0; i < delta + 1; ++i) {
    //Connects zero with third quadrant
    edges_complete.push_back({i, i + 3* nverts});
    //Connects first with second quadrant
    edges_complete.push_back({i + 2 * nverts, i + nverts});
    //Conects zero and first quadrant
    edges_complete.push_back({offset, offset + nverts});
    //Connects second and third quadrant
    edges_complete.push_back({offset + 2*nverts, offset + 3*nverts});

    offset += delta + 2 - (i+1);
  }

  side_points_complete.resize(4 * (delta+1));
  sides_complete.resize(4 * (delta+1), std::vector<bool>(4, false));
  int offset_side_pt = -1;
  int side_points_index = 0;
  //sign vector
  int index = 0;
  for (int limit = delta; limit >= 0; --limit) {
    //Assign side points
    offset_side_pt += limit+1;
    side_points_complete[side_points_index] = offset_side_pt;
    side_points_complete[side_points_index + delta+1] = offset_side_pt + nverts;
    side_points_complete[side_points_index + 2*(delta+1)] = offset_side_pt + 2*nverts;
    side_points_complete[side_points_index + 3*(delta+1)] = offset_side_pt + 3*nverts;

    sides_complete[side_points_index][0] = true;
    sides_complete[side_points_index + delta+1][1] = true;
    sides_complete[side_points_index + 2*(delta+1)][2] = true;
    sides_complete[side_points_index + 3*(delta+1)][3] = true;
    if (limit == delta) {
      sides_complete[side_points_index][3] = true;
      sides_complete[side_points_index + delta+1][0] = true;
      sides_complete[side_points_index + 2*(delta+1)][1] = true;
      sides_complete[side_points_index + 3*(delta+1)][2] = true;
    } else if (limit == 0) {
      sides_complete[side_points_index][1] = true;
      sides_complete[side_points_index][2] = true;
      sides_complete[side_points_index + delta+1][3] = true;
      sides_complete[side_points_index + 2*(delta+1)][0] = true;
    }


    side_points_index++;
    for (int y = 0; y <= limit; ++y) {
      int x = delta - limit;
      sign_complete[index] = sign_vector[index];
      sign_complete[index + nverts] = (y % 2 == 0) ? sign_vector[index] : !sign_vector[index];
      sign_complete[index + 2 * nverts] = (x % 2 == 0) ? (y % 2 == 0 ? sign_vector[index] : !sign_vector[index]) : (y % 2 == 0 ? !sign_vector[index] : sign_vector[index]);
      sign_complete[index + 3 * nverts] = (x % 2 == 0) ? sign_vector[index] : !sign_vector[index];
    ++index;
    }
  }
  // Print adjacency list for each point using edges_complete
  std::vector<std::vector<int>> adjacency_list(4 * nverts);
  for (const auto& edge : edges_complete) {
    if (sign_complete[edge.first] != sign_complete[edge.second]) { continue; }
    adjacency_list[edge.first].push_back(edge.second);
    adjacency_list[edge.second].push_back(edge.first); // assuming undirected edges
  }
}

void Graph::connected_components() {
  size_t nverts = (delta + 1) * (delta + 2) / 2;
  int n = 4 * nverts;
  std::vector<int> parent(n);
  for (int i = 0; i < n; ++i) parent[i] = i;

  auto find = [&](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };

  auto unite = [&](int x, int y) {
    int px = find(x), py = find(y);
    if (px != py) {
      parent[py] = px;
    }
  };

  // Prepare for adjacency
  std::unordered_map<int, bool> root_to_color;

  // Union vertices with same sign along the edge, otherwise record adjacency
  for (size_t i = 0; i < edges_complete.size(); ++i) {
    int u = edges_complete[i].first;
    int v = edges_complete[i].second;
    if (sign_complete[u] == sign_complete[v]) {
      unite(u, v);
    }
  }

  size_t component_count = 0;
  std::unordered_map<int, int> root_to_component;
  component.resize(n, -1);
  for (int i = 0; i < n; ++i) {
    int root = find(i);
    if (root_to_component.find(root) == root_to_component.end()) {
      root_to_component[root] = component_count++;
    }
    component[i] = root_to_component[root];
  }

  int rootptr = -1; 
  int nsidepoints = 4 * (delta+1); // total number of indices
  for (size_t i = 0; i < side_points_complete.size(); ++i) {
    int pt1 = side_points_complete[i];
    int pt2 = side_points_complete[(i+2*(delta+1)) % nsidepoints];
    int comp1 = find(pt1);
    int comp2 = find(pt2);
    if(comp1 == comp2) {
      rootptr = comp1;
      break; // Conflict found, exit loop
    }
    if (root_to_color.find(comp1) == root_to_color.end()) {
      if (root_to_color.find(comp2) == root_to_color.end()) {
        root_to_color[comp1] = true;
        root_to_color[comp2] = false;
      } else {
        root_to_color[comp1] = !root_to_color[comp2];
      }
    } else {
      if (root_to_color.find(comp2) == root_to_color.end()) {
        root_to_color[comp2] = !root_to_color[comp1];
      } else if ( root_to_color[comp2] == root_to_color[comp1]) {
      rootptr = comp1; 
      }
    }
  }

  if (rootptr >= 0) {
    root = component[rootptr];
  } else {
    std::cout << "something went wrong\n";
  }

  // Build component adjacency
  component_adjacency.clear();
  component_adjacency.resize(component_count);
  for (size_t i = 0; i < edges_complete.size(); ++i) {
    int u = edges_complete[i].first;
    int v = edges_complete[i].second;
    int cu = component[u];
    int cv = component[v];
    if (cu != cv) {
      component_adjacency[cu].insert(cv);
      component_adjacency[cv].insert(cu);
    }
  }
}

int Graph::isotopy_root() {
  connected_components();
  return root; 
}

void Graph::calculate_regions() {

  // Build component_antipod_adjacency
  std::vector<std::set<int>> component_antipode_adjacency(component_adjacency.size());
  for (size_t i = 0; i < side_points_complete.size(); ++i) {
    int comp = component[side_points_complete[i]];
    int antipodal_comp = component[side_points_complete[(i + 2 * (delta+1)) % side_points_complete.size()]];
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
        std::map<std::string, int> type_counts;
        for (const auto& child : non_leaf_children) {
          child_types.push_back(dfs(child, visited));
          type_counts[child_types.back()]++;
        }
        std::set<std::string> grouped_types;
        for (const auto& pair : type_counts) {
          grouped_types.insert(std::to_string(pair.second) + pair.first);
        }
        std::string result = open_delim;
        if (leaf_count > 0) {
          result += std::to_string(leaf_count) + sep;
        }
        //Loop of count_type_pairs to ensure order
        for (auto it = grouped_types.begin(); it != grouped_types.end(); ++it) {
          if (it != grouped_types.begin() ) { result += sep; }
          result += *it;
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

