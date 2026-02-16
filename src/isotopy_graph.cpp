#include "isotopy_graph.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cctype>
#include <mutex>

//For the Utils
#include <cstdint>

// Using declarations for common std types to reduce verbosity
using std::vector;
using std::string;
using std::pair;
using std::array;
using std::map;
using std::set;
using std::function;
using std::sort;
using std::unique;
using std::swap;
using std::make_pair;
using std::to_string;
using std::mutex;
using std::lock_guard;

// Regex types
using std::stringstream;
using std::getline;
using std::stoi;

// Exceptions
using std::out_of_range;

namespace Isotopy {

void Graph::initialize(const vector<bool>& sign_vector) {
  delta_even  = (delta % 2 == 0);
  nverts      = num_vertices(delta);
  ntotalverts = num_total_vertices(delta);
  ntriangles  = num_triangles(delta);

  // assert(sign_vector.size() == static_cast<size_t>(nverts) && "sign vector length does not match number of vertices");
  
  int q2_offset = nverts;
  int q3_offset = q2_offset + nverts - delta - 1;
  int q4_offset = q3_offset + nverts - delta - 1;

  polarisation.assign(sign_vector.begin(), sign_vector.end());
  polarisation.resize(ntotalverts, 0);
  quad_idxs.resize(nverts);
  antipodal_partner.resize(2 * (delta + 1));

  int idx = 0;
  for (int limit = delta; limit >= 0; --limit) {
    int x = delta - limit;

    for (int y = 0; y <= limit; ++y) {
      int q1_idx = idx;                                   // (+x, +y)
      int q2_idx = q2_offset + idx - delta - 1;           // (-x, +y)
      int q3_idx = q3_offset + idx - x - 1;               // (-x, -y)
      int q4_idx = q4_offset + idx - delta - 1 - x;       // (+x, -y)

      quad_idxs[idx][0] = q1_idx;

      bool q1_sign = polarisation[idx];
      bool q2_sign = (x % 2 == 0) ? q1_sign : !q1_sign;
      bool q3_sign = ((x % 2) == (y % 2)) ? q1_sign : !q1_sign;
      bool q4_sign = (y % 2 == 0) ? q1_sign : !q1_sign;

      if (x == 0 && y == 0) {
        quad_idxs[idx][1] = q1_idx;
        quad_idxs[idx][2] = q1_idx;
        quad_idxs[idx][3] = q1_idx;

      } else if (x == 0) {
        quad_idxs[idx][1] = q1_idx;
        quad_idxs[idx][2] = q3_idx;
        quad_idxs[idx][3] = q3_idx;

        polarisation[q3_idx] = q3_sign;

      } else if (y == 0) {
        quad_idxs[idx][1] = q2_idx;
        quad_idxs[idx][2] = q2_idx;
        quad_idxs[idx][3] = q1_idx;

        polarisation[q2_idx] = q2_sign;

      } else {
        quad_idxs[idx][1] = q2_idx;
        quad_idxs[idx][2] = q3_idx;
        quad_idxs[idx][3] = q4_idx;

        polarisation[q2_idx] = q2_sign;
        polarisation[q3_idx] = q3_sign;
        polarisation[q4_idx] = q4_sign;

      }

      ++idx; // increment index
    }

    // Technically two duplicate edges in antipodal_partner
    const array<int, 4>& aidxs = quad_idxs[idx - 1];
    antipodal_partner[2 * x]     = {aidxs[0], aidxs[2]};  // Q1 ↔ Q3
    antipodal_partner[2 * x + 1] = {aidxs[1], aidxs[3]};  // Q2 ↔ Q4
  }

  adjacency_edges.reserve(4 * ntriangles);
  parent.resize(ntotalverts);
  for (int i = 0; i < ntotalverts; ++i) parent[i] = i;
  rank.assign(ntotalverts, 0);
}

void Graph::process_triangles(const vector<Triangle>& triangles) {
  // assert(triangles.size() == static_cast<size_t>(ntriangles) && "Triangulation must have delta^2 triangles");

  for (const auto& [v0, v1, v2] : triangles) {
    for (int q = 0; q < 4; ++q) {
      int qv0 = quad_idxs[v0][q];
      int qv1 = quad_idxs[v1][q];
      int qv2 = quad_idxs[v2][q];

      bool s0 = polarisation[qv0];
      bool s1 = polarisation[qv1];
      bool s2 = polarisation[qv2];

      bool v0_eq_v1 = (s0 == s1);
      bool v0_eq_v2 = (s0 == s2);

      if (v0_eq_v1 && v0_eq_v2) {
        // All same sign: merge all three
        unite(parent, rank, qv0, qv1);
        unite(parent, rank, qv0, qv2);

      } else if (v0_eq_v1) {
        unite(parent, rank, qv0, qv1);
        adjacency_edges.push_back({qv0, qv2});

      } else if (v0_eq_v2) {
        unite(parent, rank, qv0, qv2);
        adjacency_edges.push_back({qv0, qv1});

      } else {
        // v1_eq_v2 (only remaining case with booleans)
        unite(parent, rank, qv1, qv2);
        adjacency_edges.push_back({qv0, qv2});
      }
    }
  }
}

Graph::Graph(int delta, const vector<bool>& sign_vector, const vector<Triangle>& triangles)
  : delta(delta) {
  initialize(sign_vector);
  process_triangles(triangles);
}

void Graph::process_edges(const vector<Edge>& edges) {
  for (const auto& [v0, v1] : edges) {
    for (int q = 0; q < 4; ++q) {
      int qv0 = quad_idxs[v0][q];
      int qv1 = quad_idxs[v1][q];

      if (polarisation[qv0] == polarisation[qv1]) {
        unite(parent, rank, qv0, qv1);
      } else {
        adjacency_edges.push_back({qv0, qv1});
      }
    }
  }
}

// Edge-based constructor
Graph::Graph(int delta, const vector<bool>& sign_vector, const vector<Edge>& edges)
  : delta(delta) {
  initialize(sign_vector);
  process_edges(edges);
}

// Backwards compatibility constructor: converts set<pair<int,int>> to vector<Edge>
Graph::Graph(int delta, const vector<bool>& sign_vector, const set<pair<int, int>>& edges)
: Graph(delta, sign_vector, vector<Edge>(edges.begin(), edges.end())) {}

// Backwards compatibility constructor: converts set<set<int>> to vector<Triangle>
Graph::Graph(int delta, const vector<bool>& sign_vector, const set<set<int>>& triangles)
: Graph(delta, sign_vector, [&triangles]() {
    vector<Triangle> tri_vec;
    tri_vec.reserve(triangles.size());
    for (const auto& tri_set : triangles) {
      // assert(tri_set.size() == 3 && "Each triangle must have exactly 3 vertices");
      vector<int> tri_tmp(tri_set.begin(), tri_set.end());
      tri_vec.push_back({tri_tmp[0], tri_tmp[1], tri_tmp[2]});
    }
    return tri_vec;
  }()) {}

void Graph::connected_components() {
  //If already computed return
  if (!component.empty()) {
    return;
  }

  // Unions already performed during process_triangles/process_edges.
  // Just assign component IDs.
  component.assign(ntotalverts, -1);
  ncomponents = 0;

  for (int i = 0; i < ntotalverts; ++i) {
    int root = find(parent, i);

    if (component[root] == -1) {
      component[root] = ncomponents++;
    }

    component[i] = component[root];
  }
}

void Graph::isotopy_type() {
  // If already done return
  if (p_regions != 0 && n_regions != -1) {
    return;
  }

  connected_components();

  // Initialize region counters
  // For odd degree, n_regions starts at 0 to account for border region at infinity
  // For even degree, n_regions starts at -1 and first increment makes it 0
  p_regions = 0;
  n_regions = delta_even ? -1 : 0;

  // Use union-find to merge antipodal components into regions
  vector<int> region_parent(ncomponents);
  for (int i = 0; i < ncomponents; ++i) region_parent[i] = i;

  // Bipartiteness check using vector indexed by component IDs (only for even degree)
  vector<int> comp_to_color;
  if (delta_even) {
    comp_to_color.assign(ncomponents, -1);  // -1 = unassigned, 0 = false, 1 = true
  }

  // Process antipodal pairs: merge regions and check bipartiteness for even degree
  if (delta_even) {
    // Even degree: merge regions AND find root via bipartiteness check
    for (const auto& [pt1, pt2] : antipodal_partner) {
      if (parent[pt1] == -1 || parent[pt2] == -1) continue;

      int root1 = find(parent, pt1);
      int root2 = find(parent, pt2);
      int comp1 = component[root1];
      int comp2 = component[root2];

      if (comp1 == comp2) {
        // Antipodal pair in same component = conflict
        if (root == -1) root = comp1;
        continue;
      }

      unite(region_parent, comp1, comp2);

      if (root == -1) {
        int& c1 = comp_to_color[comp1];
        int& c2 = comp_to_color[comp2];

        if (c1 == -1 && c2 == -1) {
          c1 = 1;
          c2 = 0;
        } else if (c1 == -1) {
          c1 = 1 - c2;
        } else if (c2 == -1) {
          c2 = 1 - c1;
        } else if (c1 == c2) {
          // Both colored with same color = conflict (odd cycle)
          root = comp1;
        }
      }
    }
  } else {
    // Odd degree: only merge regions (root found later)
    for (const auto& [pt1, pt2] : antipodal_partner) {
      if (parent[pt1] == -1 || parent[pt2] == -1) continue;

      int root1 = find(parent, pt1);
      int root2 = find(parent, pt2);
      int comp1 = component[root1];
      int comp2 = component[root2];

      if (comp1 != comp2) {
        unite(region_parent, comp1, comp2);
      }
    }
  }

  // Assign region IDs
  region.assign(ncomponents, -1);
  region_count = 0;
  for (int i = 0; i < ncomponents; ++i) {
    int region_root = find(region_parent, i);
    if (region[region_root] == -1) {
      region[region_root] = region_count++;
    }
    region[i] = region[region_root];
  }

  // Build region adjacency directly from adjacency_edges
  region_adjacency.clear();
  region_adjacency.resize(region_count);

  if (!delta_even && root == -1) {
    // Odd degree with no root yet: find self-loop in region adjacency
    for (const auto& [u, v] : adjacency_edges) {
      int cu = component[u];
      int cv = component[v];
      if (cu != cv) {
        int rcu = region[cu];
        int rcv = region[cv];
        region_adjacency[rcu].push_back(rcv);
        region_adjacency[rcv].push_back(rcu);

        if (rcu == rcv) {
          root = cu;
        }
      }
    }
  } else {
    // Even degree or root already found: just build adjacency
    for (const auto& [u, v] : adjacency_edges) {
      int cu = component[u];
      int cv = component[v];
      if (cu != cv) {
        int rcu = region[cu];
        int rcv = region[cv];
        region_adjacency[rcu].push_back(rcv);
        region_adjacency[rcv].push_back(rcu);
      }
    }
  }

  for (auto& nbrs : region_adjacency) {
    sort(nbrs.begin(), nbrs.end());
    nbrs.erase(unique(nbrs.begin(), nbrs.end()), nbrs.end());
  }

  if (!delta_even) {
    // Increment region count by 1 for the border region
    region_count++;
  }

  root_region = region[root];

  // BFS to assign region signs (2-coloring of region adjacency graph)
  region_sign.assign(region_count, false);
  region_sign[root_region] = !delta_even;  // true for odd degree, false for even

  vector<bool> visited(region_count, false);
  visited[root_region] = true;

  vector<int> queue;
  queue.reserve(region_count);
  queue.push_back(root_region);

  for (size_t i = 0; i < queue.size(); ++i) {
    int current = queue[i];

    if (region_sign[current]) {
      p_regions++;
    } else {
      n_regions++;
    }

    for (int neighbor : region_adjacency[current]) {
      if (!visited[neighbor]) {
        visited[neighbor] = true;
        region_sign[neighbor] = !region_sign[current];
        queue.push_back(neighbor);
      }
    }
  }

}

string Graph::viro_notation(bool unicode) {
  if (region_adjacency.empty()) {
    isotopy_type();
  }

  //Im odd delta fall sollte bei der fake region angefangen werden
  if (!delta_even) {
    string notation = Isotopy::viro_notation(region[root], region_adjacency, unicode);
    const string open_delim = unicode ? "\u27E8" : "<";
    const string close_delim = unicode ? "\u27E9" : ">";
    const bool has_additional_components = notation.size() > open_delim.size() + close_delim.size();
    const string insert_fragment = has_additional_components ? (unicode ? "J\u2294" : "Jv") : "J";
    notation.insert(open_delim.size(), insert_fragment);
    return notation;
  }
  return Isotopy::viro_notation(region[root], region_adjacency, unicode);
}

string viro_notation(int root_region, const Adjacency& region_adjacency, bool unicode) {
  const string open_delim = unicode ? "\u27E8" : "<";
  const string close_delim = unicode ? "\u27E9" : ">";
  const string sep = unicode ? "\u2294" : "v";

  vector<bool> visited(region_adjacency.size(), false);

  function<string(int)> dfs =
    [&](int curr_region) -> string {
      visited[curr_region] = true;
      int leaf_count = 0;
      vector<int> non_leaf_children;
      for (const auto& neighbor : region_adjacency[curr_region]) {
        if (!visited[neighbor]) {
          bool is_leaf = true;
          for (const auto& nn : region_adjacency[neighbor]) {
            if (!visited[nn] && nn != curr_region) {
              is_leaf = false;
              break;
            }
          }
          if (is_leaf) {
            leaf_count++;
            visited[neighbor] = true;
          } else {
            non_leaf_children.push_back(neighbor);
          }
        }
      }
      if (non_leaf_children.empty()) {
        if (leaf_count == 0) {
          return open_delim  + close_delim;
        }
        return open_delim + to_string(leaf_count) + close_delim;
      } else if (leaf_count == 0 && non_leaf_children.size() == 1) {
        return open_delim + "1" + dfs(non_leaf_children[0]) + close_delim;
      } else {
        vector<string> child_types;
        map<string, int> type_counts;
        for (const auto& child : non_leaf_children) {
          child_types.push_back(dfs(child));
          type_counts[child_types.back()]++;
        }
        vector<pair<int,string>> grouped_types;
        for (const auto& pair : type_counts) {
          grouped_types.push_back(make_pair(pair.second, pair.first));
        }
        sort(grouped_types.begin(), grouped_types.end(), [](const pair<int, string>& a, const pair<int, string>& b) {
            if (a.second.length() != b.second.length()) {
            return a.second.length() < b.second.length();
            }
            return a.second < b.second;
            });

        string result = open_delim;
        if (leaf_count > 0) {
          result += to_string(leaf_count) + sep;
        }
        //Loop of count_type_pairs to ensure order
        for (auto it = grouped_types.begin(); it != grouped_types.end(); ++it) {
          if (it != grouped_types.begin()) { result += sep; }
          result += to_string(it->first) + it->second;
        }
        result += close_delim;
        return result;
      }
    };

  return dfs(root_region);
}

int num_vertices(int delta) {
  return (delta + 1) * (delta + 2) / 2;
}

int num_total_vertices(int delta) {
  return 2 * delta * delta + 2 * delta + 1;
}

int num_edges(int delta) {
  return 3 * delta + (3 * (delta * delta - delta)) / 2;
}

int num_triangles(int delta) {
  return delta * delta;
}

vector<Edge> triangles_to_edges(const vector<Triangle>& triangles) {
  vector<Edge> edges;
  edges.reserve(triangles.size() * 3);
  for (const auto& tri : triangles) {
    const int v0 = tri[0];
    const int v1 = tri[1];
    const int v2 = tri[2];
    auto push_edge = [&](int a, int b) {
      if (a > b) swap(a, b);
      edges.emplace_back(a, b);
    };
    push_edge(v0, v1);
    push_edge(v0, v2);
    push_edge(v1, v2);
  }
  sort(edges.begin(), edges.end());
  edges.erase(unique(edges.begin(), edges.end()), edges.end());
  return edges;
}

vector<Edge> triangles_to_edges(const std::set<std::set<int>>& triangles) {
  vector<Triangle> tri_vec;
  tri_vec.reserve(triangles.size());
  for (const auto& tri_set : triangles) {
    vector<int> sorted_tri(tri_set.begin(), tri_set.end());
    if (sorted_tri.size() != 3) {
      throw std::invalid_argument("Each triangle must have exactly 3 vertices.");
    }
    tri_vec.push_back({sorted_tri[0], sorted_tri[1], sorted_tri[2]});
  }
  return triangles_to_edges(tri_vec);
}

// Helper to build coordinate mapping for a given delta (cached)
// This duplicates the quadrant indexing logic from Graph::initialize
// but provides a standalone version for point_to_idx/idx_to_point functions
static const vector<pair<int,int>>& get_coord_table(int delta) {
  static map<int, vector<pair<int,int>>> cache;
  static mutex cache_mutex;
  lock_guard<mutex> lock(cache_mutex);

  auto it = cache.find(delta);
  if (it != cache.end()) {
    return it->second;
  }

  int nverts = num_vertices(delta);
  int ntotal = num_total_vertices(delta);
  int q2_offset = nverts;
  int q3_offset = q2_offset + nverts - delta - 1;
  int q4_offset = q3_offset + nverts - delta - 1;

  vector<pair<int,int>> coords(ntotal);

  int idx = 0;
  for (int limit = delta; limit >= 0; --limit) {
    int x = delta - limit;
    for (int y = 0; y <= limit; ++y) {
      coords[idx] = {x, y};

      if (x == 0 && y == 0) {
        // Origin: all quadrants map to same point
      } else if (x == 0) {
        // Y-axis: Q1/Q2 same, Q3/Q4 same
        int q3_idx = q3_offset + idx - x - 1;
        coords[q3_idx] = {-x, -y};
      } else if (y == 0) {
        // X-axis: Q1/Q4 same, Q2/Q3 same
        int q2_idx = q2_offset + idx - delta - 1;
        coords[q2_idx] = {-x, y};
      } else {
        // General case: all four quadrants distinct
        int q2_idx = q2_offset + idx - delta - 1;
        int q3_idx = q3_offset + idx - x - 1;
        int q4_idx = q4_offset + idx - delta - 1 - x;
        coords[q2_idx] = {-x, y};
        coords[q3_idx] = {-x, -y};
        coords[q4_idx] = {x, -y};
      }
      ++idx;
    }
  }

  cache[delta] = std::move(coords);
  return cache[delta];
}

int point_to_idx(int delta, int x, int y) {
  // Build reverse lookup cache: (x,y) -> idx
  static map<int, map<pair<int,int>, int>> reverse_cache;
  static mutex cache_mutex;
  lock_guard<mutex> lock(cache_mutex);

  auto delta_it = reverse_cache.find(delta);
  if (delta_it == reverse_cache.end()) {
    // Build reverse mapping from coord table
    const auto& coords = get_coord_table(delta);
    auto& reverse_map = reverse_cache[delta];
    for (size_t i = 0; i < coords.size(); ++i) {
      reverse_map[coords[i]] = i;
    }
    delta_it = reverse_cache.find(delta);
  }

  auto it = delta_it->second.find({x, y});
  return (it != delta_it->second.end()) ? it->second : -1;
}

pair<int,int> idx_to_point(int delta, int idx) {
  int ntotal = num_total_vertices(delta);
  if (idx < 0 || idx >= ntotal) {
    throw out_of_range("idx_to_point: index out of range");
  }
  const auto& coords = get_coord_table(delta);
  return coords[idx];
}

}  // namespace Isotopy

namespace Utils {

// map<pair<int,int>, int> get_pt2int(int delta) {
//   map<pair<int,int>, int> pt2int;
//   int index = 0;
//   size_t nverts = Isotopy::num_vertices(delta);
//   for (int orthant = 0; orthant < 4; ++orthant) {
//     for (int limit = delta; limit >= 0; --limit) {
//       for (int y = 0; y <= limit; ++y) {
//         int x = delta - limit;
//         int adjusted_y = (orthant == 1 || orthant == 2) ? -y : y;
//         int adjusted_x = (orthant == 2 || orthant == 3) ? -x : x;
//         pt2int[{adjusted_x, adjusted_y}] = index + orthant * nverts;
//         ++index;
//       }
//     }
//     index = 0; // Reset index for the next orthant
//   }

//   return pt2int;
// }

pair<vector<bool>, vector<Isotopy::Triangle>> pcom_to_signs_and_triangles_vec(const string& pcom_string_input) {
  vector<bool> sign_vector;
  vector<Isotopy::Triangle> triangles;

  string pcom_string = pcom_string_input;
  const string cols_token = ",{\"cols\":";
  size_t cols_pos = 0;
  while ((cols_pos = pcom_string.find(cols_token, cols_pos)) != string::npos) {
    size_t end = pcom_string.find('}', cols_pos);
    if (end == string::npos) break;
    pcom_string.erase(cols_pos, end - cols_pos + 1);
  }

  auto bracket_range = [&](const string& key) -> pair<size_t, size_t> {
    size_t start = pcom_string.find(key);
    if (start == string::npos) throw std::runtime_error("Missing key in pcom string: " + key);
    size_t open = pcom_string.find('[', start);
    if (open == string::npos) throw std::runtime_error("Malformed pcom: '[' not found for " + key);
    int depth = 0;
    size_t pos = open;
    for (; pos < pcom_string.size(); ++pos) {
      if (pcom_string[pos] == '[') ++depth;
      else if (pcom_string[pos] == ']') --depth;
      if (depth == 0) break;
    }
    if (depth != 0) throw std::runtime_error("Malformed pcom: unmatched brackets for " + key);
    return {open, pos};
  };

  auto [sign_open, sign_close] = bracket_range("\"SIGNS\"");
  for (size_t i = sign_open + 1; i < sign_close; ) {
    char c = pcom_string[i];
    if (c == 't') {
      sign_vector.push_back(true);
      i += 4;
    } else if (c == 'f') {
      sign_vector.push_back(false);
      i += 5;
    } else {
      ++i;
    }
  }

  auto [cells_open, cells_close] = bracket_range("\"MAXIMAL_CELLS\"");
  size_t cursor = cells_open + 1;
  while (true) {
    cursor = pcom_string.find('[', cursor);
    if (cursor == string::npos || cursor >= cells_close) break;
    ++cursor;

    Isotopy::Triangle tri{};
    for (int idx = 0; idx < 3; ++idx) {
      while (cursor < cells_close &&
             !(pcom_string[cursor] == '-' || isdigit(static_cast<unsigned char>(pcom_string[cursor]))))
        ++cursor;
      size_t end = cursor;
      while (end < cells_close &&
             (pcom_string[end] == '-' || isdigit(static_cast<unsigned char>(pcom_string[end]))))
        ++end;
      if (cursor == end) throw std::runtime_error("Malformed triangle data in pcom");
      tri[idx] = stoi(pcom_string.substr(cursor, end - cursor));
      cursor = pcom_string.find_first_of(",]", end);
      if (cursor == string::npos || cursor > cells_close)
        throw std::runtime_error("Malformed triangle entry in pcom");
      if (pcom_string[cursor] == ',') ++cursor;
    }
    triangles.push_back(tri);
    cursor = pcom_string.find(']', cursor);
    if (cursor == string::npos || cursor > cells_close)
      throw std::runtime_error("Malformed triangle closing bracket in pcom");
    ++cursor;
  }

  return make_pair(sign_vector, triangles);
}

std::pair<vector<bool>, std::set<std::set<int>>> pcom_to_signs_and_triangles(const string& pcom_string) {
  auto [signs, triangles_vec] = pcom_to_signs_and_triangles_vec(pcom_string);
  std::set<std::set<int>> triangles_set;
  for (const auto& tri : triangles_vec) {
    triangles_set.insert({tri[0], tri[1], tri[2]});
  }
  return {signs, triangles_set};
}

string signs_and_triangles_to_pcom(const vector<bool>& sign_vector, const vector<Isotopy::Triangle>& triangles) {
    return signs_and_triangles_to_pcom(sign_vector, triangles, "<unknown>");
}
string signs_and_triangles_to_pcom(const vector<bool>& sign_vector, const vector<Isotopy::Triangle>& triangles, string origin_tag) {
  const int lib_vers = 2; 
  string faces_str = "\"MAXIMAL_CELLS\": ";
  faces_str += "[\n";
  for (const auto& triangle : triangles) {
    //faces_str += triangulation.dcel.to_string(face) + ",\n";
    faces_str += "[";
    for (auto point : triangle) {
        faces_str += to_string(point);
        faces_str += ","; 
    }
    faces_str.pop_back();
    faces_str += "],\n";
  }
  faces_str.pop_back();
  faces_str.pop_back();
  faces_str += "]\n";
  string sign_str = "\"SIGNS\": [\n";
  for (size_t i = 0; i < sign_vector.size(); ++i) {
    string sign = sign_vector[i] ? "true" : "false";
    sign_str += sign + ",\n";
  }
  sign_str.pop_back();
  sign_str.pop_back();
  sign_str += "]\n";
  string filet_str = "{ \"_ns\": { \"polymake\": [ \"https://polymake.org\", \"4.13\" ] },\n";
  filet_str += "  \"_type\": \"tropical::Hypersurface<Min>\",\n";
  filet_str += "  \"_id\": \"filet\",\n";
  filet_str += "  \"_libisotopy_version\": \""+ to_string(lib_vers) +"\",\n"; 
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

std::pair<vector<bool>, std::set<std::set<int>>> pcom_to_signs_and_triangles_set(const string& pcom_string) {
  return pcom_to_signs_and_triangles(pcom_string);
}

string signs_and_triangles_to_pcom(const vector<bool>& sign_vector, const std::set<std::set<int>>& triangles) {
  return signs_and_triangles_to_pcom(sign_vector, triangles, "<unknown>");
}

string signs_and_triangles_to_pcom(const vector<bool>& sign_vector, const std::set<std::set<int>>& triangles, string origin_tag) {
  // Convert set<set<int>> to vector<Triangle>
  vector<Isotopy::Triangle> triangles_vec;
  triangles_vec.reserve(triangles.size());
  for (const auto& tri_set : triangles) {
    vector<int> tri_tmp(tri_set.begin(), tri_set.end());
    triangles_vec.push_back({tri_tmp[0], tri_tmp[1], tri_tmp[2]});
  }
  return signs_and_triangles_to_pcom(sign_vector, triangles_vec, origin_tag);
}

std::map<std::pair<int,int>, int> get_pt2int(int delta) {
  std::map<std::pair<int,int>, int> pt2int;
  int nverts = Isotopy::num_total_vertices(delta);
  for (int idx = 0; idx < nverts; ++idx) {
    auto [x, y] = Isotopy::idx_to_point(delta, idx);
    pt2int[{x, y}] = idx;
  }
  return pt2int;
}

}  // namespace Utils
