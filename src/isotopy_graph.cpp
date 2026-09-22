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
#include <numeric>

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
  if (delta > MAX_DELTA) {
    throw std::invalid_argument("delta=" + std::to_string(delta) + " exceeds MAX_DELTA=" + std::to_string(MAX_DELTA));
  }

  delta_even  = (delta % 2 == 0);
  nverts      = num_vertices(delta);
  ntotalverts = num_total_vertices(delta);
  ntriangles  = num_triangles(delta);

  if (sign_vector.size() != static_cast<size_t>(nverts)) {
    throw std::invalid_argument("sign vector has " + std::to_string(sign_vector.size()) +
                                " entries, expected " + std::to_string(nverts) +
                                " for delta=" + std::to_string(delta));
  }
  
  int q2_offset = nverts;
  int q3_offset = q2_offset + nverts - delta - 1;
  int q4_offset = q3_offset + nverts - delta - 1;

  for (int i = 0; i < nverts; ++i) polarisation[i] = sign_vector[i];  // Copy Q1 signs (Q2-Q4 set below by quadrant reflection)
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
  for (int i = 0; i < ntotalverts; ++i) parent[i] = i;  // Initialize union-find: each vertex is its own parent
  std::fill(rank.begin(), rank.begin() + ntotalverts, 0);  // Zero-initialize rank for used portion only
}

void Graph::process_triangles(const vector<Triangle>& triangles) {
  if (triangles.size() != static_cast<size_t>(ntriangles)) {
    throw std::invalid_argument("triangulation has " + std::to_string(triangles.size()) +
                                " triangles, expected " + std::to_string(ntriangles) +
                                " for delta=" + std::to_string(delta));
  }

  for (const auto& [v0, v1, v2] : triangles) {
    assert(v0 >= 0 && v0 < nverts && v1 >= 0 && v1 < nverts && v2 >= 0 && v2 < nverts &&
           "triangle vertex index out of range");
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
    assert(v0 >= 0 && v0 < nverts && v1 >= 0 && v1 < nverts &&
           "edge vertex index out of range");
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
      if (tri_set.size() != 3) {
        throw std::invalid_argument("Each triangle must have exactly 3 vertices.");
      }
      vector<int> tri_tmp(tri_set.begin(), tri_set.end());
      tri_vec.push_back({tri_tmp[0], tri_tmp[1], tri_tmp[2]});
    }
    return tri_vec;
  }()) {}

void Graph::connected_components() {
  //If already computed return
  if (components_computed) {
    return;
  }

  // Unions already performed during process_triangles/process_edges.
  // Just assign component IDs.
  std::fill(component.begin(), component.begin() + ntotalverts, -1);  // Fill only used portion of fixed-size array
  ncomponents = 0;

  for (int i = 0; i < ntotalverts; ++i) {
    int root = find(parent, i);

    if (component[root] == -1) {
      component[root] = ncomponents++;
    }

    component[i] = component[root];
  }
  components_computed = true;
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
  array<int, MAX_TOTAL_VERTS> region_parent;  // Stack-allocated union-find for regions
  for (int i = 0; i < ncomponents; ++i) region_parent[i] = i;

  // Bipartiteness check using array indexed by component IDs (only for even degree)
  array<int, MAX_TOTAL_VERTS> comp_to_color;  // Stack-allocated color array
  if (delta_even) {
    std::fill(comp_to_color.begin(), comp_to_color.begin() + ncomponents, -1);  // -1 = unassigned, 0 = false, 1 = true
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
  std::fill(region.begin(), region.begin() + ncomponents, -1);  // Fill only used portion of fixed-size array
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

  // Two distinct components can share a region, so rcu == rcv is possible. Such an
  // edge is not an adjacency between regions; recording it would put a self-loop in
  // region_adjacency, which breaks the leaf test in viro_notation(). For odd degree
  // it is instead the signal that identifies the isotopy root.
  if (!delta_even && root == -1) {
    // Odd degree with no root yet: find self-loop in region adjacency
    for (const auto& [u, v] : adjacency_edges) {
      int cu = component[u];
      int cv = component[v];
      if (cu != cv) {
        int rcu = region[cu];
        int rcv = region[cv];
        if (rcu == rcv) {
          root = cu;
        } else {
          region_adjacency[rcu].push_back(rcv);
          region_adjacency[rcv].push_back(rcu);
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
        if (rcu != rcv) {
          region_adjacency[rcu].push_back(rcv);
          region_adjacency[rcv].push_back(rcu);
        }
      }
    }
  }

  for (auto& nbrs : region_adjacency) {
    sort(nbrs.begin(), nbrs.end());
    nbrs.erase(unique(nbrs.begin(), nbrs.end()), nbrs.end());
  }

  if (!delta_even) {
    // Increment region count by 1 for the border region, and give it an (empty)
    // adjacency slot so that region_adjacency.size() == region_count always holds.
    // The border region is unreachable from root_region, so the sign BFS below and
    // the p_regions / n_regions tallies are unaffected; n_regions is seeded to 0
    // for odd degree to account for it.
    region_count++;
    region_adjacency.resize(region_count);
  }

  if (root == -1) {
    throw std::runtime_error("no isotopy root found for delta=" + std::to_string(delta) +
                             "; the patchwork is not a valid T-curve");
  }

  root_region = region[root];

  // BFS to assign region signs (2-coloring of region adjacency graph)
  std::fill(region_sign.begin(), region_sign.begin() + region_count, 0);  // Zero-initialize only used portion
  region_sign[root_region] = !delta_even;  // true for odd degree, false for even

  array<uint8_t, MAX_TOTAL_VERTS> visited;  // Stack-allocated visited flags
  std::fill(visited.begin(), visited.begin() + region_count, 0);  // Zero-initialize only used portion
  visited[root_region] = 1;

  array<int, MAX_TOTAL_VERTS> bfs_queue;  // Stack-allocated BFS queue (replaces vector with push_back)
  int bfs_size = 0;
  bfs_queue[bfs_size++] = root_region;

  for (int i = 0; i < bfs_size; ++i) {
    int current = bfs_queue[i];

    if (region_sign[current]) {
      p_regions++;
    } else {
      n_regions++;
    }

    for (int neighbor : region_adjacency[current]) {
      if (!visited[neighbor]) {
        visited[neighbor] = 1;
        region_sign[neighbor] = !region_sign[current];
        bfs_queue[bfs_size++] = neighbor;
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
            // nn != neighbor: a caller-supplied adjacency may carry self-loops, which
            // are not children and must not defeat the leaf test.
            if (!visited[nn] && nn != curr_region && nn != neighbor) {
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

vector<RegionTreeNode> get_region_tree(int root_region, const Adjacency& region_adjacency) {
  vector<RegionTreeNode> result;
  vector<bool> visited(region_adjacency.size(), false);

  function<void(int, int)> dfs = [&](int curr, int parent_idx) {
    int my_idx = (int)result.size();
    result.push_back({curr, parent_idx, 0});
    visited[curr] = true;

    int leaf_count = 0;
    vector<int> non_leaf_children;
    for (int nb : region_adjacency[curr]) {
      if (!visited[nb]) {
        bool is_leaf = true;
        for (int nn : region_adjacency[nb])
          if (!visited[nn] && nn != curr && nn != nb) { is_leaf = false; break; }
        if (is_leaf) { leaf_count++; visited[nb] = true; }
        else           non_leaf_children.push_back(nb);
      }
    }
    result[my_idx].leaf_count = leaf_count;
    for (int child : non_leaf_children) dfs(child, my_idx);
  };

  dfs(root_region, -1);
  return result;
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

// Reconstruct triangulation faces from an edge list.
// Finds all 3-cliques in the adjacency graph, but only keeps unit-area triangles
// (|area2|==1 by Pick's theorem) to avoid false positives where three vertices are
// pairwise connected but a fourth vertex sits inside their convex hull.
vector<Triangle> edges_to_triangles(const vector<Edge>& edges, int delta) {
  map<int, set<int>> adj;
  for (auto [u, v] : edges) { adj[u].insert(v); adj[v].insert(u); }

  set<Triangle> seen;
  vector<Triangle> tris;
  for (auto [u, v] : edges) {
    for (int w : adj[u]) {
      if (adj[v].count(w)) {
        Triangle t = {u, v, w};
        sort(t.begin(), t.end());
        if (seen.insert(t).second) {
          auto [xa, ya] = idx_to_point(delta, t[0]);
          auto [xb, yb] = idx_to_point(delta, t[1]);
          auto [xc, yc] = idx_to_point(delta, t[2]);
          int area2 = (xb-xa)*(yc-ya) - (yb-ya)*(xc-xa);
          if (area2 == 1 || area2 == -1)
            tris.push_back(t);
        }
      }
    }
  }
  return tris;
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

int original_idx(int delta, int idx) {
  const auto& coords = get_coord_table(delta);
  auto [x, y] = coords[idx];
  return point_to_idx(delta, abs(x), abs(y));
}

void Graph::apply_flip(const Edge& old_edge, const Edge& new_edge) {
  auto replace_or_move_edge = [](vector<pair<int,int>>& from, vector<pair<int,int>>& to, Edge old_edge, Edge new_edge, bool replace_only) {
    auto is_old_edge = [&](const pair<int,int>& e) {
      return (e == old_edge) || (e.first == old_edge.second && e.second == old_edge.first);
    };

    bool found = false;
    for (auto& e : from) {
      if (is_old_edge(e)) {
        found = true;
        if (replace_only) {
          e = new_edge;
        }
      }
    }
    if (!replace_only) {
      from.erase(std::remove_if(from.begin(), from.end(), is_old_edge), from.end());
      to.push_back(new_edge);
    } else if (!found) {
      from.push_back(new_edge);
    }
  };

  for (int q = 0; q < 4; ++q) {
    int old_v0 = quad_idxs[old_edge.first][q];
    int old_v1 = quad_idxs[old_edge.second][q];
    bool old_s0 = polarisation[old_v0];
    bool old_s1 = polarisation[old_v1];
    int new_v0 = quad_idxs[new_edge.first][q];
    int new_v1 = quad_idxs[new_edge.second][q];
    bool new_s0 = polarisation[new_v0];
    bool new_s1 = polarisation[new_v1];

    if (old_s0 == old_s1 && new_s0 == new_s1) {
      // Both are component edges: replace
      replace_or_move_edge(component_edges, adjacency_edges, Edge{old_v0, old_v1}, Edge{new_v0, new_v1}, true);
    } else if (old_s0 != old_s1 && new_s0 != new_s1) {
      // Both are adjacency edges: replace
      replace_or_move_edge(adjacency_edges, component_edges, Edge{old_v0, old_v1}, Edge{new_v0, new_v1}, true);
    } else if (old_s0 == old_s1 && new_s0 != new_s1) {
      // Old is component, new is adjacency: move
      replace_or_move_edge(component_edges, adjacency_edges, Edge{old_v0, old_v1}, Edge{new_v0, new_v1}, false);
    } else if (old_s0 != old_s1 && new_s0 == new_s1) {
      // Old is adjacency, new is component: move
      replace_or_move_edge(adjacency_edges, component_edges, Edge{old_v0, old_v1}, Edge{new_v0, new_v1}, false);
    }
  }
  
  invalidate_cache();
}

void Graph::invalidate_cache() {
  // Reset union-find and rebuild from current triangulation + signs.
  // master's connected_components() relies on unions pre-built by process_triangles,
  // so we must re-run it here rather than leaving parent in identity state.
  component_edges.clear();
  adjacency_edges.clear();
  for (int i = 0; i < ntotalverts; ++i) { parent[i] = i; rank[i] = 0; }
  if (!tri_list.empty()) process_triangles(tri_list);
  // Reset result caches
  components_computed = false;
  region_adjacency.clear();
  p_regions = 0;
  n_regions = -1;
  root = -1;
  root_region = -1;
  region_count = 0;
}

void Graph::prepare(const vector<Triangle>& triangles) {

  tri_list = triangles;

  // Triangles need to be counter clockwise oriented
  for (auto& tri : tri_list) {
    auto [x0, y0] = idx_to_point(delta, tri[0]);
    auto [x1, y1] = idx_to_point(delta, tri[1]);
    auto [x2, y2] = idx_to_point(delta, tri[2]);
    int cross = (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
    if (cross < 0) swap(tri[1], tri[2]);
  }

  int n = static_cast<int>(tri_list.size());
  half_edges.resize(3 * n);
  vertex_he.assign(nverts, -1);

  // Initialise each triangle's three half-edges
  for (int i = 0; i < n; ++i) {
    auto [v0, v1, v2] = tri_list[i];
    half_edges[3*i]   = {v0, -1, 3*i+1, 3*i+2, i};
    half_edges[3*i+1] = {v1, -1, 3*i+2, 3*i+0, i};
    half_edges[3*i+2] = {v2, -1, 3*i+0, 3*i+1, i};
  }

  he_map.clear();
  for (int k = 0; k < 3*n; ++k) {
    int dest = half_edges[half_edges[k].next].origin;
    he_map[{half_edges[k].origin, dest}] = k;
  }

  // Link twins
  for (int k = 0; k < 3*n; ++k) {
    int dest = half_edges[half_edges[k].next].origin;
    auto it = he_map.find({dest, half_edges[k].origin});
    if (it != he_map.end()) {
      half_edges[k].twin = it->second;
    }
  }

  // One outgoing half-edge per vertex
  for (int k = 0; k < 3*n; ++k) {
    vertex_he[half_edges[k].origin] = k;
  }

  // Build normalized edge set
  edge_set.clear();
  for (int k = 0; k < 3*n; ++k) {
    int u = half_edges[k].origin;
    int v = half_edges[half_edges[k].next].origin;
    edge_set.insert(u < v ? Edge{u, v} : Edge{v, u});
  }

  // Verify: all different-sign (adjacency) edges must appear in the supplied triangulation.
  for (const auto& e : adjacency_edges) {
    int u = original_idx(delta, e.first);
    int v = original_idx(delta, e.second);
    Edge norm = (u < v) ? Edge{u, v} : Edge{v, u};
    if (!edge_set.count(norm))
      throw std::invalid_argument("prepare(): triangulation edge not found in current graph structure");
  }

  // I cleanly rebuild the component_edges and adjacency_edges sets, pretty sure that this is unnecessary could be removed in the future.
  set<Edge> new_component_edges;
  set<Edge> new_adjacency_edges;
  for (const auto& tri : tri_list) {
    for (auto [a, b] : {pair{tri[0],tri[1]}, pair{tri[1],tri[2]}, pair{tri[2],tri[0]}}) {
      for (int q = 0; q < 4; ++q) {
        int qa = quad_idxs[a][q];
        int qb = quad_idxs[b][q];
        int sign_a = polarisation[qa];
        int sign_b = polarisation[qb];
        Edge norm_q = (qa < qb) ? Edge{qa, qb} : Edge{qb, qa};
      if (sign_a == sign_b) {
        new_component_edges.insert(norm_q);
      } else {
        new_adjacency_edges.insert(norm_q);
      }
      }
    }
  }
  component_edges.assign(new_component_edges.begin(), new_component_edges.end());
  adjacency_edges.assign(new_adjacency_edges.begin(), new_adjacency_edges.end());

}

bool Graph::is_flippable(const Edge& edge) {
  // Find correct half-edge
  auto it = he_map.find({edge.first, edge.second});
  if (it == he_map.end()) it = he_map.find({edge.second, edge.first});
  if (it == he_map.end()) return false;

  int hei = it->second;
  int tw  = half_edges[hei].twin;
  if (tw < 0) return false;
  if (half_edges[hei].face < 0 || half_edges[tw].face < 0) return false;

  // The 4 vertices of the quad, ordered [u, op1, v, op2] for convexity check
  int u   = half_edges[hei].origin;
  int v   = half_edges[tw].origin;
  int op1 = half_edges[half_edges[hei].prev].origin;
  int op2 = half_edges[half_edges[tw].prev].origin;

  auto [x1, y1] = idx_to_point(delta, u);
  auto [x2, y2] = idx_to_point(delta, op1);
  auto [x3, y3] = idx_to_point(delta, v);
  auto [x4, y4] = idx_to_point(delta, op2);
  array<pair<int,int>, 4> pts = {{{x1,y1},{x2,y2},{x3,y3},{x4,y4}}};

  // All four points must be distinct
  for (int i = 0; i < 4; ++i)
    for (int j = i+1; j < 4; ++j)
      if (pts[i] == pts[j]) return false;

  // No collinear triple, and all cross products must have the same sign (convex quad)
  int prev_cross = 0;
  for (int i = 0; i < 4; ++i) {
    int j = (i+1)%4, k = (i+2)%4;
    int area = pts[i].first * (pts[j].second - pts[k].second)
             + pts[j].first * (pts[k].second - pts[i].second)
             + pts[k].first * (pts[i].second - pts[j].second);
    if (area == 0) return false;
    if (prev_cross == 0) prev_cross = area;
    else if ((area > 0) != (prev_cross > 0)) return false;
  }
  return true;
}

void Graph::update_dcel(int hei, int tw, int u, int v, int w, int x) {
  // Named indices for the 6 half-edges involved in the flip:
  //   hei = u→v (being replaced by x→w)
  //   tw  = v→u (being replaced by w→x)
  //   bc  = v→w (next of hei)
  //   ca  = w→u (prev of hei)
  //   ad  = u→x (next of tw)
  //   db  = x→v (prev of tw)
  int bc = half_edges[hei].next;
  int ca = half_edges[hei].prev;
  int ad = half_edges[tw].next;
  int db = half_edges[tw].prev;

  int f1 = half_edges[hei].face;
  int f2 = half_edges[tw].face;

  // Repurpose hei as x→w and tw as w→x
  half_edges[hei].origin = x;
  half_edges[hei].next   = ca;
  half_edges[hei].prev   = ad;
  // twin and face stay the same

  half_edges[tw].origin = w;
  half_edges[tw].next   = db;
  half_edges[tw].prev   = bc;
  // twin and face stay the same

  // Relink the four surrounding half-edges
  half_edges[bc].next = tw;  half_edges[bc].prev = db;
  half_edges[ca].next = ad;  half_edges[ca].prev = hei;
  half_edges[ad].next = hei; half_edges[ad].prev = ca;
  half_edges[db].next = bc;  half_edges[db].prev = tw;

  // ad and bc change faces
  half_edges[ad].face = f1;
  half_edges[bc].face = f2;

  // Update tri_list
  tri_list[f1] = {u, x, w};
  tri_list[f2] = {v, w, x};

  // Update vertex_he: u and v lost their outgoing half-edge
  if (vertex_he[u] == hei) vertex_he[u] = ad;
  if (vertex_he[v] == tw)  vertex_he[v] = bc;
  vertex_he[x] = hei;  // hei now has origin x
  vertex_he[w] = tw;   // tw now has origin w

  // Update he_map
  he_map.erase({u, v});
  he_map.erase({v, u});
  he_map[{x, w}] = hei;
  he_map[{w, x}] = tw;

  // Update edge_set
  Edge old_norm = (u < v) ? Edge{u, v} : Edge{v, u};
  Edge new_norm = (w < x) ? Edge{w, x} : Edge{x, w};
  edge_set.erase(old_norm);
  edge_set.insert(new_norm);
}

vector<Edge> Graph::flippable_edges() {
  vector<Edge> result;
  for (const auto& e : edge_set)
    if (is_flippable(e)) result.push_back(e);
  return result;
}

void Graph::update_edge(const Edge& edge) {
  if (!is_flippable(edge)) return;

  auto it = he_map.find({edge.first, edge.second});
  if (it == he_map.end()) it = he_map.find({edge.second, edge.first});

  int hei = it->second;
  int tw  = half_edges[hei].twin;
  int u   = half_edges[hei].origin;
  int v   = half_edges[tw].origin;
  int w   = half_edges[half_edges[hei].prev].origin;
  int x   = half_edges[half_edges[tw].prev].origin;

  update_dcel(hei, tw, u, v, w, x);   
  invalidate_cache();                  
}

void Graph::update_sign(int vector_index) {
  // Why this exists: Boundary vertices can share the same global index across quadrants
  std::set<int> flipped;
  for (int q = 0; q < 4; ++q) {
    int qv = quad_idxs[vector_index][q];
    if (flipped.insert(qv).second)
      polarisation[qv] = !polarisation[qv];
  }
  invalidate_cache();
}

}  // namespace Isotopy

namespace Utils {

vector<int> transpose_partition(const vector<int>& input){
   int size= input.back();
   vector<int> result(size);
   for(const auto& i: input){
      for(int j=0; j<i; j++){
         result[size-1-j]++;
      }
   }
   return result;
}

vector<vector<int>> partitions_of_height_k(int n, int k){
   vector<vector<int>> result;
   if(n<k){ return result; }
   if(k==1){
      vector<int> tmp(n);
      for(unsigned int i=0; i<tmp.size(); i++){
         tmp[i]=1;
      }
      result.push_back(tmp);
      return result;
   }
   if(k == n){
      vector<int> tmp{n};
      result.push_back(tmp);
      return result;
   }
   for(int i=1; i<=k; i++){
      auto tail = partitions_of_height_k(n-k, i);
      for(const auto& p: tail){
         vector<int> tmp(p);
         tmp.push_back(k);
         result.push_back(tmp);
      }
   }
   return result;
}

vector<vector<int>> partitions_of_length_k(int n, int k){
   vector<vector<int>> result;
   for(const auto& p : partitions_of_height_k(n, k)){
      result.push_back(transpose_partition(p));
   }
   return result;
}

vector<vector<int>> partitions_of_max_length_k(int n, int k){
   vector<vector<int>> result;
   for(int i=1; i<=k; i++){
      auto next = partitions_of_length_k(n, i);
      result.insert(result.end(),
            std::make_move_iterator(next.begin()),
            std::make_move_iterator(next.end()));
   }
   return result;
}

vector<std::string> trees_of_size(int n){
   std::string root("o");
   if(n == 1){ return vector<std::string>{root}; }
   vector<std::string> result;
   std::string trivial = "o[";
   bool first = true;
   for(int i=1; i<=n-1; i++){
      trivial += (first ? "o" : ",o");
      first = false;
   }
   trivial += "]";
   result.push_back(trivial);
   for(int i=1; i<=n-1; i++){
      for(const auto& p: partitions_of_max_length_k(n-i-1, i)){
         auto it = p.rbegin();
         vector<std::string> inner_res{"o["};
         bool first = true;
         while(it != p.rend()){
            auto sts = trees_of_size(*it+1);
            vector<std::string> newres{};
            for(const auto& t : sts){
               for(const auto& r : inner_res){
                  std::string concat = r + (first ? "" : ",") + t;
                  newres.push_back(concat);
               }
            }
            std::swap(newres, inner_res);
            ++it;
            first = false;
         }
         std::string tail = "";
         vector<std::string> newres;
         first = true;
         for(int j=0; j<i-p.size(); j++){
            tail += (first ? "o" : ",o");
            first = false;
         }
         tail = tail + "]";
         for(const auto& r : inner_res){
            std::string concat = r + (i>p.size() ? "," : "")+ tail;
            newres.push_back(concat);
         }
         std::swap(newres, inner_res);
         result.insert(result.end(),
               std::make_move_iterator(inner_res.begin()),
               std::make_move_iterator(inner_res.end()));
      }
   }
   return result;
}

// Copied from hilbert_encoding.pl
int orientation(int i, int j, int k, const map<int, vector<int>>& int2pt){
  const vector<int>& p0(int2pt.at(i)), p1(int2pt.at(j)), p2(int2pt.at(k));
  int det=0;
  for(int a = 0; a<3; a++){
    det += p0[a] * p1[(a+1)%3] * p2[(a+2)%3];
    det -= p0[a] * p1[(a+2)%3] * p2[(a+1)%3];
  }
  if(det == 0) return 0;
  else if(det >0) return 1;
  else return -1;
}

bool segments_intersect(int i, int j, int k, int l, const map<int, vector<int>>& int2pt){
  int ijk = orientation(i, j, k, int2pt);
  if(ijk == 0) return false;
  int ijl = orientation(i, j, l, int2pt);
  if(ijk == 0 || ijk == ijl) return false;
  int kli = orientation(k, l, i, int2pt);
  if(kli ==0) return false;
  int klj = orientation(k, l, j, int2pt);
  if(klj == 0 || kli == klj) return false;
  return true;
}


void remove_intersecting_edges(vector<pair<int, int>>& edges, const pair<int, int>& edge, const map<int, vector<int>>& int2pt){
  edges.erase(
    std::remove_if(edges.begin(), edges.end(),
        [&](const std::pair<int,int>& x) {
            return segments_intersect(x.first, x.second, edge.first, edge.second, int2pt);
        }),
    edges.end()
);
}

pair<int, int> get_and_remove_random_edge(vector<pair<int, int>>& edges){
  pair<int, int> result;
  int pos = std::rand() % edges.size();
  auto iter = edges.begin();
  for(int i=0; i<pos; i++){
    iter++;
  }
  result = *iter;
  edges.erase(iter);
  return result;
}

vector<pair<int, int>> get_random_triangulation(int delta, const vector<pair<int, int>>& existing_edges) {
// void get_random_triangulation(int delta) {
  map<vector<int>, int> pt2int;
  map<int, vector<int>> int2pt;
  int ptct = 0;
  for(int x = 0; x<=delta; x++){
    for(int y = 0; y<=delta-x; y++){
      vector<int> pt{x,y,1};
      pt2int[pt] = ptct;
      int2pt[ptct] = pt;
      ptct++;
    }
  }
  vector<pair<int, int>> primitive_edges;
  for(int i=0; i<ptct; i++){
    for(int j=i+1; j<ptct; j++){
      vector<int> p0(int2pt[i]), p1(int2pt[j]);
      int check = std::gcd(p0[0]-p1[0], p0[1]-p1[1]);
      if(check == 1){
        primitive_edges.push_back(pair<int, int>{i, j});
      }
    }
  }
  for(const auto& e: existing_edges){
    remove_intersecting_edges(primitive_edges, e, int2pt);
  }
  // std::cout << "Length after init: " << primitive_edges.size() << std::endl;
  vector<pair<int, int>> result;
  while(!primitive_edges.empty()){
    pair<int, int> ne(get_and_remove_random_edge(primitive_edges));
    result.push_back(ne);
    // std::cout << "Adding edge: " << ne.first << " " << ne.second << std::endl;
    // std::cout << "Result length: " << result.size() << std::endl;
    remove_intersecting_edges(primitive_edges, ne, int2pt);
  }
  // int target = 3 * delta + (3 * (delta * delta - delta)) / 2;
  // std::cout << "Target is " << target << std::endl;
  return result;
}

vector<pair<int, int>> get_random_triangulation(int delta) {
  vector<pair<int, int>> ee;
  return get_random_triangulation(delta, ee);
}

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
