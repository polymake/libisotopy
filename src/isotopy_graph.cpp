#include "isotopy_graph.h"
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <functional>
#include <regex>
#include <string>
#include <sstream>

//For the Utils
#include <cstdint>

namespace Isotopy {

Graph::Graph(int delta, const std::vector<bool>& sign_vector, const std::set<std::set<int>>& triangles) {
  size_t nverts = (delta + 1) * (delta + 2) / 2;
  size_t nedges = 3 * delta + (3 * (delta * delta - delta)) / 2;

  size_t upper_size = nverts * (nverts - 1) / 2;
  std::vector<bool> adjacency(upper_size, false);

  auto idx = [nverts](size_t i, size_t j) {
    if (i > j) std::swap(i, j);
    return i * nverts - (i * (i + 1)) / 2 + (j - i - 1);
  };

  for (const auto& triangle : triangles) {
    if (triangle.size() != 3) {
      throw std::invalid_argument("Each triangle must have exactly 3 vertices.");
    }
    auto it1 = triangle.begin();
    for (size_t a = 0; a < 2; ++a, ++it1) {
      auto it2 = std::next(it1);
      for (size_t b = a + 1; b < 3; ++b, ++it2) {
        adjacency[idx(*it1, *it2)] = true;
      }
    }
  }

  std::vector<std::pair<int, int>> edges;
  edges.reserve(nedges);
  for (size_t i = 0; i < nverts; ++i) {
    for (size_t j = i + 1; j < nverts; ++j) {
      if (adjacency[idx(i, j)]) {
        edges.push_back({static_cast<int>(i), static_cast<int>(j)});
      }
    }
  }
  *this = Graph(delta, sign_vector, edges);
}

Graph::Graph(int delta, const std::vector<bool>& sign_vector, const std::set<std::pair<int, int>>& edges)
: Graph(delta, sign_vector, std::vector<std::pair<int, int>>(edges.begin(), edges.end())) {}


Graph::Graph(int delta, const std::vector<bool>&sign_vector, const std::vector<std::pair<int, int>>& edges)
: delta(delta) {
  size_t nverts = (delta + 1) * (delta + 2) / 2;
  assert(sign_vector.size() == nverts && "sign vector length does not match number of vertices");

  sign_complete.resize(4 * nverts);
  edges_complete.reserve(4 * edges.size() + 4*(delta+1));
  for (size_t i = 0; i < 4; ++i) {
    for (const auto& edge : edges) {
      edges_complete.push_back({edge.first + i * nverts, edge.second + i * nverts});
    }
  }

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
  //If already computed return
  if (!component.empty()) {
    return;
  }

  size_t nverts = (delta + 1) * (delta + 2) / 2;
  int n = 4 * nverts;
  parent.resize(n);
  for (int i = 0; i < n; ++i) parent[i] = i;


  // Union vertices with same sign along the edge, otherwise record adjacency
  for (size_t i = 0; i < edges_complete.size(); ++i) {
    int u = edges_complete[i].first;
    int v = edges_complete[i].second;
    if (sign_complete[u] == sign_complete[v]) {
      unite(u, v);
    }
  }

  //This is completely useless except for the fact, that it counts the components 
  size_t component_count = 0;
  std::unordered_map<int, int> root_to_component;
  component.resize(n, -1);
  for (int i = 0; i < n; ++i) {
    int component_root_ptr = find(i);
    if (root_to_component.find(component_root_ptr) == root_to_component.end()) {
      root_to_component[component_root_ptr] = component_count++;
    }
    component[i] = root_to_component[component_root_ptr];
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
  //if already done return
  if (root != -1) {
    return root;
  }
  if (component.empty()) {
    connected_components();
  }

  if (delta % 2 == 0) {

    // Prepare for adjacency
    std::unordered_map<int, bool> root_to_color;

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
    }
   } else {
    calculate_regions();

    //increment region + 1 since, border is one
    region_count++;

    root_region = -1;
    for (size_t i = 0; i < component_adjacency.size(); ++i) {
      for (int adj_comp : component_adjacency[i]) {
        if (region[i] == region[adj_comp]) {
          root_region = region[i];
          root = i;
          break;
        }
      }
      if (root_region != -1) break;
    }
  }
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
  region_count = 0;
  for (size_t i = 0; i < region.size(); ++i) {
    if (region[i] == -1) {
      std::set<int> stack;
      stack.insert(i);
      while (!stack.empty()) {
        int current = *stack.begin();
        stack.erase(stack.begin());
        region[current] = region_count;
        for (int neighbor : component_antipode_adjacency[current]) {
          if (region[neighbor] == -1) {
            stack.insert(neighbor);
          }
        }
      }
      region_count++;
    }
  }


  region_adjacency.resize(region_count, std::set<int>());
  for (size_t i = 0; i < component_adjacency.size(); ++i) {
    for (int adj_comp : component_adjacency[i]) {
      region_adjacency[region[i]].insert(region[adj_comp]);
    }
  }

}

void Graph::isotopy_type() {
  //If already done return
  if (p_regions != 0 && n_regions != -1) {
    return;
  }

  if (root == -1) {
    isotopy_root();
  }

  if (delta % 2 == 0) {
    calculate_regions();
  }

  //Might cause issues in the future

  root_region = region[root];

  region_sign.resize(region_adjacency.size(), false);

  std::set<std::pair<int, int>> edges;

  std::set<int> stack;
  stack.insert(root_region);
  if (delta % 2 != 0) {
    region_sign[root_region] = true;
    n_regions++;
  } else {
    region_sign[root_region] = false; 
  }

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

  //Im odd delta fall sollte bei der fake region angefangen werden
  if (delta % 2 == 1) {
    std::string notation = Isotopy::viro_notation(region[root], region_adjacency, unicode);
    const std::string open_delim = unicode ? "\u27E8" : "<";
    const std::string close_delim = unicode ? "\u27E9" : ">";
    const bool has_additional_components = notation.size() > open_delim.size() + close_delim.size();
    const std::string insert_fragment = has_additional_components ? (unicode ? "J\u2294" : "Jv") : "J";
    notation.insert(open_delim.size(), insert_fragment);
    return notation;
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
        if (leaf_count == 0) {
          return open_delim  + close_delim;
        }
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
        std::vector<std::pair<int,std::string>> grouped_types;
        for (const auto& pair : type_counts) {
          grouped_types.push_back(std::make_pair(pair.second, pair.first));
        }
        std::sort(grouped_types.begin(), grouped_types.end(), [](const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
            if (a.second.length() != b.second.length()) {
            return a.second.length() < b.second.length();
            }
            return a.second < b.second;
            });

        std::string result = open_delim;
        if (leaf_count > 0) {
          result += std::to_string(leaf_count) + sep;
        }
        //Loop of count_type_pairs to ensure order
        for (auto it = grouped_types.begin(); it != grouped_types.end(); ++it) {
          if (it != grouped_types.begin()) { result += sep; }
          result += std::to_string(it->first) + it->second;
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

std::map<std::pair<int,int>, int> get_pt2int(int delta) {
  std::map<std::pair<int,int>, int> pt2int;
  int index = 0;
  int nverts = (delta + 1) * (delta + 2) / 2;
  for (int orthant = 0; orthant < 4; ++orthant) {
    for (int limit = delta; limit >= 0; --limit) {
      for (int y = 0; y <= limit; ++y) {
        int x = delta - limit;
        int adjusted_y = (orthant == 1 || orthant == 2) ? -y : y;
        int adjusted_x = (orthant == 2 || orthant == 3) ? -x : x;
        pt2int[{adjusted_x, adjusted_y}] = index + orthant * nverts;
        ++index;
      }
    }
    index = 0; // Reset index for the next orthant
  }

  return pt2int;
}


std::pair<std::vector<bool>, std::set<std::set<int>>> pcom_to_signs_and_triangles(const std::string& pcom_string_input) {
    std::vector<bool> sign_vector;
    std::set<std::set<int>> triangles;

    //bring files in common format
    std::regex cols_regex(",\\{\"cols\":[0-9]+\\}");
    std::string cols_deleted = std::regex_replace(pcom_string_input, cols_regex, "");
    std::regex whitespaces_regex("\\s");
    std::string pcom_string = std::regex_replace(cols_deleted, whitespaces_regex, "");
    //std::cout << pcom_string << std::endl; 

    //extract sign vector
    std::regex signs_regex("\"SIGNS\":\\[([^ \\]]*)\\]");
    std::smatch base_match; 
    std::regex_search(pcom_string, base_match, signs_regex);
    std::stringstream base_match_stream(base_match[1]);
    std::string temp_sign;  
    char del = ','; 
    while (getline(base_match_stream, temp_sign, del)){
            if (temp_sign == "true"){
                sign_vector.push_back(1);
            } else {
                    sign_vector.push_back(0);
            }
    }
    /*for (auto i : sign_vector){
        std::cout << i;
    }
    std::cout << std::endl; */
        

    //extract triangles
    std::regex triangle_regex("\\[([0-9]+,[0-9]+,[0-9]+)\\]");
    std::smatch triangle_match;
    std::regex_search(pcom_string, triangle_match, triangle_regex);
    auto triangles_begin = std::sregex_iterator(pcom_string.begin(), pcom_string.end(), triangle_regex);
    auto triangles_end = std::sregex_iterator();
    for(std::sregex_iterator i = triangles_begin; i!=triangles_end; ++i){
        std::smatch match = *i;
        std::string temp_tri = match[1].str(); 
        std::stringstream temp_tri_stream(temp_tri); 
        std::string temp_number; 
        char del = ','; 
        std::set<int> triangle; 
        while (getline(temp_tri_stream, temp_number, del)){
                triangle.insert(std::stoi(temp_number));
        }
        //std::cerr << match.str() << std::endl;
        triangles.insert(triangle);
    }
/*
    for (auto i : triangles){
        for (auto x:i){
            std::cout << x << ",";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl; */
    

    return std::make_pair(sign_vector, triangles);
}

std::string signs_and_triangles_to_pcom(const std::vector<bool>& sign_vector, const std::set<std::set<int>>& triangles) {
    return signs_and_triangles_to_pcom(sign_vector, triangles, "<unknown>");
}
std::string signs_and_triangles_to_pcom(const std::vector<bool>& sign_vector, const std::set<std::set<int>>& triangles, std::string origin_tag) {
  const int lib_vers = 2; 
  std::string faces_str = "\"MAXIMAL_CELLS\": ";
  faces_str += "[\n";
  for (const auto& triangle : triangles) {
    //faces_str += triangulation.dcel.to_string(face) + ",\n";
    faces_str += "[";
    for (auto point : triangle) {
        faces_str += std::to_string(point) += ","; 
    }
    faces_str.pop_back();
    faces_str += "],\n";
  }
  faces_str.pop_back();
  faces_str.pop_back();
  faces_str += "]\n";
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
  filet_str += "  \"_libisotopy_version\": \""+ std::to_string(lib_vers) +"\",\n"; 
  filet_str += "  \"_attrs\": { \"ORIGIN\": { \"attachment\": true } },\n";
  filet_str += "  \"ORIGIN\": \"" + origin_tag + "\",\n";
  filet_str += "  \"DUAL_SUBDIVISION\": {\n";
  filet_str += faces_str;
  filet_str += "},\n";
  filet_str += "\"PATCHWORK\": [{\n";
  filet_str += "  \"_id\": \"filet#0\",\n";
  filet_str += sign_str;
  filet_str += "}]}\n";
  return filet_str;
}

}
