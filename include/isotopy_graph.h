#pragma once

#include <vector>
#include <set>
#include <map>
#include <string>

namespace Isotopy {

/**
 * @brief Represents an isotopy type as a graph and provides methods for isotopy type computation.
 *
 * The Graph struct encapsulates the data and algorithms needed to construct a graph from
 * side points and edges or triangles, compute isotopy invariants, and generate Viro notation.
 */
struct Graph {
  int delta; ///< The degree of the patchwork.
  std::vector<bool> sign; ///< Sign vector with a boolean value for each vertex.
  std::vector<std::set<int>> adjacency; ///< Adjacency list for the graph.

  std::vector<int> side_points; ///< Maps side point index to vertex index.
  std::vector<std::pair<int,int>> int2pt; ///< Maps vertex index to (x,y) coordinates in the triangle.

  // Information about connected components
  int root = -1; ///< Root component index.
  std::vector<int> component; ///< component[i] gives the component index of vertex i.
  std::vector<std::set<int>> component_adjacency; ///< component_adjacency[c] gives the set of components adjacent to component c.

  int root_region = -1; ///< Region index of the root component.
  int region_count = 0; ///< Total number of regions.
  std::vector<int> region; ///< region[c] gives the region index of component c.
  std::vector<std::set<int>> region_adjacency; ///< region_adjacency[r] gives the set of regions adjacent to region r.
  std::vector<bool> region_sign; ///< region_sign[r] gives the sign of region r (true for positive, false for negative).
  

  std::set<std::pair<int, int>> edges = std::set<std::pair<int, int>>(); ///< Edges of the tree representing the isotopy graph.

  std::vector<int> side_points_complete;
  std::vector<std::pair<int,int>> edges_complete;
  std::vector<bool> sign_complete;
  std::vector<std::vector<bool>> sides_complete;

  int p_regions = 0; ///< Number of even regions.
  int n_regions = -1; ///< Number of odd regions.


  /**
   * @brief Constructs an empty graph.
   *
   * Initializes an empty graph with no vertices, edges, or triangles.
   */
  Graph() = default;

  /**
   * @brief Constructs a graph from a sign vector and a set of edges.
   *
   * @param delta The degree of the patchwork.
   * @param sign_vector The sign vector with a boolean value for each vertex.
   * @param edges The set of edges, each represented as a pair of vertex indices.
   */
  Graph(int delta, const std::vector<bool>& sign_vector, const std::set<std::pair<int, int>>& edges);


  /**
   * @brief Constructs a graph from a sign vector and a list of edges.
   *
   * @param delta The degree of the patchwork.
   * @param sign_vector The sign vector with a boolean value for each vertex.
   * @param edges The list of edges, each represented as a pair of vertex indices.
   */
  Graph(int delta, const std::vector<bool>& sign_vector, const std::vector<std::pair<int, int>>& edges);



  /**
   * @brief Constructs a graph from a sign vector and a set of triangles.
   *
   * @param delta defines the degree of the patchwork
   * @param sign_vector The sign vector with a boolean value for each vertex.
   * @param triangles The set of triangles, each represented as a set of three vertex indices.
   */
  Graph(int delta, const std::vector<bool>& sign_vector, const std::set<std::set<int>>& triangles);

  /**
   * @brief Computes the connected components of the graph.
   *
   * Populates the component-related member variables, assigning each vertex to a component.
   * Updates the component adjacency information.
   */
  void connected_components();

  /**
   * @brief Prepares neighbor sets for isotopy root computation.
   *
   * For each side point, computes the set of other side points it is connected to via its component.
   * Populates the `nbs` member with this neighbor information.
   */
  void pre_isotopy_root();

  /**
   * @brief Computes the isotopy root of the graph.
   *
   * This is the core function of the library. It iteratively merges side points and their antipodes
   * according to the isotopy rules, updating neighbor and component information until a root is found
   * or a maximum number of iterations is reached.
   *
   * The isotopy root is defined as the component containing a side point that is connected to its antipode,
   * indicating that the isotopy process has merged all relevant regions.
   *
   * @return The component index of the isotopy root if found, or -1 if no root is found after the maximum iterations.
   */
  int isotopy_root();

  /**
   * @brief Calculates the regions of the graph based on component antipode adjacency.
   *
   * Assigns each component to a region by traversing antipodal connections, grouping components
   * that are connected via antipodes into the same region. Also builds the region adjacency structure
   * and determines the root region.
   *
   * This function must be called after the isotopy root has been computed.
   */
  void calculate_regions();

  /**
   * @brief Computes the full isotopy type of the graph.
   *
   * This is the recommended entry point for users. It performs all necessary steps to analyze the graph's
   * isotopy type, including computing connected components, finding the isotopy root, calculating regions,
   * and assigning region signs. After calling this function, all relevant isotopy invariants and structures
   * are available for further queries.
   *
   * @throws std::runtime_error if no isotopy root is found.
   */
  void isotopy_type();

  /**
   * @brief Returns the Viro notation for the isotopy type of the graph.
   *
   * Generates a string representation (Viro notation) of the isotopy type, using the computed regions and their
   * adjacencies. Optionally outputs Unicode symbols for delimiters and separators.
   *
   * @param unicode If true, use Unicode symbols; otherwise, use ASCII.
   * @return The Viro notation string for the isotopy type.
   */
  std::string viro_notation(bool unicode = false);

  /**
   * @brief Returns the number of even regions in the isotopy decomposition.
   *
   * @return The count of even regions.
   */
  int even_regions() {
    if (root_region == -1) isotopy_type();
    return p_regions; };

  /**
   * @brief Returns the number of odd regions in the isotopy decomposition.
   *
   * @return The count of odd regions.
   */
  int odd_regions() {
    if (root_region == -1) isotopy_type(); 
    return n_regions; 
  };

  /**
   * @brief Checks if the graph represents an M-curve.
   *
   * An M-curve is characterized by having the maximal possible number of regions for the given degree.
   *
   * @return True if the graph is an M-curve, false otherwise.
   */
  bool is_mcurve() {
    if (root_region == -1) isotopy_type();
    int expected_regions = (delta - 1) * (delta - 2) / 2 + 1; 
    return (p_regions + n_regions) == expected_regions;
  };

  void lazy_compute(bool need_components, bool need_root, bool need_regions, bool need_region_counts);

  std::vector<int> parent;

  int find(int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  }

  void unite(int x, int y) {
    int px = find(x), py = find(y);
    if (px != py) {
      parent[py] = px;
    }
  }

};

/**
 * @brief Generates the Viro notation string for a given root region and region adjacency structure.
 *
 * This function can be used independently to produce Viro notation from explicit region and adjacency data.
 *
 * @param root_region The index of the root region.
 * @param region_adjacency The adjacency list of regions.
 * @param unicode If true, use Unicode symbols; otherwise, use ASCII.
 * @return The Viro notation string.
 */
std::string viro_notation(int root_region, const std::vector<std::set<int>>& region_adjacency, bool unicode = false);

/**
 * @brief Computes the number of vertices in a triangular grid of degree delta.
 *
 * For a patchwork of degree delta, the triangular grid has (delta+1)(delta+2)/2 vertices.
 * This is the number of lattice points in a triangle with side length delta+1.
 *
 * @param delta The degree of the patchwork.
 * @return The number of vertices in the triangular grid.
 */
int num_vertices(int delta);

/**
 * @brief Computes the number of edges in a triangular grid of degree delta.
 *
 * The total edges equal the boundary edges (3*delta) plus the interior edges
 * (three per interior lattice segment), giving 3*delta + 3*(delta^2 - delta)/2.
 *
 * @param delta The degree of the patchwork.
 * @return The number of edges in the triangular grid.
 */
int num_edges(int delta);

/**
 * @brief Converts a set of triangles to a vector of edges.
 *
 * Each triangle is represented as a set of three vertex indices. This function extracts all edges
 * from the triangles, ensuring consistent ordering (smaller index first) and removing duplicates.
 * Uses an efficient adjacency matrix approach for optimal performance.
 *
 * @param delta The degree of the patchwork (determines number of vertices).
 * @param triangles A set of triangles, where each triangle is a set of three vertex indices.
 * @return A vector of edges, where each edge is a pair of vertex indices with the smaller index first.
 */
std::vector<std::pair<int, int>> triangles_to_edges(int delta, const std::set<std::set<int>>& triangles);

}

namespace Utils {

std::map<std::pair<int,int>, int> get_pt2int(int delta);

std::pair<std::vector<bool>, std::set<std::set<int>>> pcom_to_signs_and_triangles(const std::string& pcom_string);

std::string signs_and_triangles_to_pcom(const std::vector<bool>& sign_vector, const std::set<std::set<int>>& triangles);
std::string signs_and_triangles_to_pcom(const std::vector<bool>& sign_vector, const std::set<std::set<int>>& triangles, std::string origin_tag);

}
