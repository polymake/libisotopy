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

Graph::Graph(int delta, const std::vector<bool>& sign_vector, const std::set<std::pair<int, int>>& edges)
: Graph(delta, sign_vector, std::vector<std::pair<int, int>>(edges.begin(), edges.end())) {}


Graph::Graph(int delta, const std::vector<bool>&sign_vector, const std::vector<std::pair<int, int>>& edges)
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
/*
std::string create_polymake_string(int delta, const std::vector<bool>& sign_vector, const std::set<std::set<int>>& triangles) {
  auto& triangulation = this->triangulation;
  std::string faces_str = "\"MAXIMAL_CELLS\": ";
  faces_str += "[";
  for (const auto& face : triangulation.dcel.faces) {
    faces_str += triangulation.dcel.to_string(face) + ",\n";
  }
  faces_str.pop_back();
  faces_str.pop_back();
  faces_str += "]\n";
  std::string points_str = "\"POINTS\": [\n";
  for (const auto& pt : triangulation.dcel.vertices) {
    std::string fcoord = std::to_string(triangulation.delta- pt.first - pt.second);
    points_str += "[\"1\", \"" + std::to_string(pt.first) + "\", \"" + std::to_string(pt.second)+"\", \""+fcoord + "\"]" + ",\n";
  }
  points_str.pop_back();
  points_str.pop_back();
  points_str += "]\n";
  std::string sign_str = "\"SIGNS\": [\n";
  for (size_t i = 0; i < sign_vector.size(); ++i) {
    std::string sign = sign_vector[i] ? "true" : "false";
    sign_str += sign + ",\n";
  }
  sign_str.pop_back();
  sign_str.pop_back();
  sign_str += "]\n";
  std::string filet_str = "{ \"_ns\": { \"polymake\": [ \"https://polymake.org\", \"4.13\" ] },\n";
  filet_str += "  \"_type\": \"tropical::Hypersurface<Min>\",\n";
  filet_str += "  \"_id\": \"filet\",\n";
  filet_str += "  \"DUAL_SUBDIVISION\": {\n";
  filet_str += points_str + ",\n";
  filet_str += faces_str;
  filet_str += "},\n";
  filet_str += "\"PATCHWORK\": [{\n";
  filet_str += "  \"_id\": \"filet#0\",\n";
  filet_str += sign_str + "\n";
  filet_str += "  }]}\n";
  return filet_str;
}

void Patchworking::write_polymake_file(const std::string& filePath) {
  std::string json = polymake_filet();
  std::ofstream out(filePath);
  out << json;
  out.close();
}
*/
}

