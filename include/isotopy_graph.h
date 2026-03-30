#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Isotopy {

// Bring common std types into Isotopy namespace for cleaner declarations
using std::vector;
using std::string;
using std::pair;
using std::array;

// Domain-specific type aliases for clarity
using Triangle = array<int, 3>;
using Edge = pair<int, int>;
using Adjacency = vector<vector<int>>;
using QuadrantIndices = array<int, 4>;

// Maximum supported degree for fixed-size array optimization
static constexpr int MAX_VERTS = (MAX_DELTA + 1) * (MAX_DELTA + 2) / 2;
static constexpr int MAX_TOTAL_VERTS = 2 * MAX_DELTA * MAX_DELTA + 2 * MAX_DELTA + 1;

/**
 * @brief Represents an isotopy type as a graph and provides methods for isotopy type computation.
 *
 * The Graph struct encapsulates the data and algorithms needed to construct a graph from
 * side points and edges or triangles, compute isotopy invariants, and generate Viro notation.
 */
struct Graph {
  int delta = 0;        ///< The degree of the patchwork.
  bool delta_even;      ///< True if delta is even, false if odd.
  int nverts = 0;       ///< Number of vertices in the first quadrant.
  int ntotalverts = 0;  ///< Total number of compressed vertices across quadrants.
  int ntriangles = 0;   ///< Number of triangles in the input triangulation.
  int ncomponents = 0;  ///< Number of connected components (after connected_components()).

  // Quadrant-based representation
  array<uint8_t, MAX_TOTAL_VERTS> polarisation; ///< polarisation[vertex_idx] gives the sign of vertex vertex_idx (initialized in initialize())
  array<QuadrantIndices, MAX_VERTS> quad_idxs;  ///< quad_idxs[i][q] gives the global index of vertex i reflected to quadrant q (initialized in initialize())

  // Edge lists for component analysis
  vector<Edge> component_edges; ///< Edges for connected component analysis (same-sign edges, no longer populated)
  vector<Edge> adjacency_edges; ///< Edges for component adjacency (different-sign edges)

  // Connected component information
  array<int, MAX_TOTAL_VERTS> parent; ///< Union-find parent array for component construction (initialized in initialize())
  array<int, MAX_TOTAL_VERTS> rank;   ///< Union-find rank array for union-by-rank optimization (initialized in initialize())
  int root = -1; ///< Root component index.
  array<int, MAX_TOTAL_VERTS> component; ///< component[i] gives the component index of vertex i. (initialized in connected_components())
  bool components_computed = false; ///< Guard for connected_components() idempotency (replaces component.empty() check)
  Adjacency component_adjacency; ///< component_adjacency[c] gives the set of components adjacent to component c.

  // Region information
  int root_region = -1; ///< Region index of the root component.
  int region_count = 0; ///< Total number of regions.
  array<int, MAX_TOTAL_VERTS> region; ///< region[c] gives the region index of component c. (initialized in isotopy_type())
  Adjacency region_adjacency; ///< region_adjacency[r] gives the set of regions adjacent to region r.
  array<uint8_t, MAX_TOTAL_VERTS> region_sign; ///< region_sign[r] gives the sign of region r (initialized in isotopy_type(), 0=negative, 1=positive)

  int p_regions = 0; ///< Number of positive regions.
  int n_regions = -1; ///< Number of negative regions.

  // Legacy members (may be deprecated)
  vector<Edge> antipodal_partner;
  vector<Edge> edges_complete;


  /**
   * @brief Constructs an empty graph.
   *
   * Initializes an empty graph with no vertices, edges, or triangles.
   */
  Graph() = default;

  /**
   * @brief Constructs a graph with only delta set (for two-phase construction).
   *
   * Use this followed by initialize() and process_triangles()/process_edges()
   * for fine-grained control or benchmarking individual phases.
   *
   * @param delta The degree of the patchwork.
   */
  explicit Graph(int delta) : delta(delta) {}

  /**
   * @brief Initializes the graph structure from a sign vector.
   *
   * Sets up quadrant indices, polarisation array, antipodal partners, and
   * reserves space for edge vectors. Must be called before process_triangles()
   * or process_edges().
   *
   * @param sign_vector The sign vector with a boolean value for each vertex.
   */
  void initialize(const vector<bool>& sign_vector);

  /**
   * @brief Processes triangles to populate component and adjacency edges.
   *
   * Must be called after initialize(). Classifies triangle edges as same-sign
   * (component edges) or different-sign (adjacency edges) across all quadrants.
   *
   * @param triangles The list of triangles, each represented as three vertex indices.
   */
  void process_triangles(const vector<Triangle>& triangles);

  /**
   * @brief Processes edges to populate component and adjacency edges.
   *
   * Must be called after initialize(). Classifies edges as same-sign
   * (component edges) or different-sign (adjacency edges) across all quadrants.
   *
   * @param edges The list of edges, each represented as a pair of vertex indices.
   */
  void process_edges(const vector<Edge>& edges);

  /**
   * @brief Constructs a graph from a sign vector and a list of edges.
   *
   * @param delta The degree of the patchwork.
   * @param sign_vector The sign vector with a boolean value for each vertex.
   * @param edges The list of edges, each represented as a pair of vertex indices.
   */
  Graph(int delta, const vector<bool>& sign_vector, const vector<Edge>& edges);

  /**
   * @brief Constructs a graph from a sign vector and a set of edges (backwards compatibility).
   *
   * @param delta The degree of the patchwork.
   * @param sign_vector The sign vector with a boolean value for each vertex.
   * @param edges The set of edges, each represented as a pair of vertex indices.
   */
  Graph(int delta, const vector<bool>& sign_vector, const std::set<std::pair<int, int>>& edges);

  /**
   * @brief Constructs a graph from a sign vector and a list of triangles.
   *
   * @param delta defines the degree of the patchwork
   * @param sign_vector The sign vector with a boolean value for each vertex.
   * @param triangles The list of triangles, each represented as three vertex indices.
   */
  Graph(int delta, const vector<bool>& sign_vector, const vector<Triangle>& triangles);

  /**
   * @brief Constructs a graph from a sign vector and a set of triangles (backwards compatibility).
   *
   * @param delta defines the degree of the patchwork
   * @param sign_vector The sign vector with a boolean value for each vertex.
   * @param triangles The set of triangles, each represented as a set of three vertex indices.
   */
  Graph(int delta, const vector<bool>& sign_vector, const std::set<std::set<int>>& triangles);

  /**
   * @brief Computes the connected components of the graph.
   *
   * Populates the component-related member variables, assigning each vertex to a component.
   * Updates the component adjacency information.
   */
  void connected_components();

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
  string viro_notation(bool unicode = false);

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

private:
  // Union-find helpers templated to work with both array and vector
  template<typename T>
  static int find(T& parent_array, int x) {
    while (parent_array[x] != x) {
      parent_array[x] = parent_array[parent_array[x]];
      x = parent_array[x];
    }
    return x;
  }

  template<typename T>
  static void unite(T& parent_array, int x, int y) {
    int px = find(parent_array, x);
    int py = find(parent_array, y);
    if (px != py) parent_array[py] = px;
  }

  template<typename T, typename R>
  static void unite(T& parent_array, R& rank_array, int x, int y) {
    int px = find(parent_array, x);
    int py = find(parent_array, y);
    if (px != py) {
      if (rank_array[px] < rank_array[py]) {
        parent_array[px] = py;
      } else if (rank_array[px] > rank_array[py]) {
        parent_array[py] = px;
      } else {
        parent_array[py] = px;
        rank_array[px]++;
      }
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
string viro_notation(int root_region, const Adjacency& region_adjacency, bool unicode = false);

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
 * @brief Computes the total number of compressed vertices across all quadrants.
 *
 * This counts all distinct lattice points obtained by reflecting the triangular grid
 * of degree delta into the four quadrants and identifying coincident points.
 *
 * @param delta The degree of the patchwork.
 * @return The total number of distinct vertices.
 */
int num_total_vertices(int delta);

/**
 * @brief Returns the compressed vertex index for a given lattice point.
 *
 * The index matches the internal indexing used by Graph (via quad_idxs and signs),
 * obtained by reflecting the degree-delta triangular grid into all four quadrants
 * and identifying coincident points.
 *
 * @param delta The degree of the patchwork.
 * @param x The x-coordinate of the lattice point.
 * @param y The y-coordinate of the lattice point.
 * @return The vertex index, or -1 if the point is not part of the grid.
 */
int point_to_idx(int delta, int x, int y);

/**
 * @brief Returns the lattice point corresponding to a compressed vertex index.
 *
 * @param delta The degree of the patchwork.
 * @param idx The vertex index.
 * @return The (x,y) coordinates of the lattice point.
 * @throws std::out_of_range if idx is outside the valid range.
 */
pair<int,int> idx_to_point(int delta, int idx);

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
 * @brief Computes the number of triangles in a triangulation of degree delta.
 *
 * For a complete triangulation of a triangular grid with degree delta, there are delta^2 triangles.
 *
 * @param delta The degree of the patchwork.
 * @return The number of triangles in the triangulation.
 */
int num_triangles(int delta);

/**
 * @brief Converts triangles to a unique edge list.
 *
 * @param triangles Triangles represented as Isotopy::Triangle entries.
 * @return A deduplicated vector of edges.
 */
vector<Edge> triangles_to_edges(const vector<Triangle>& triangles);

/**
 * @brief Converts set-based triangles to a unique edge list.
 *
 * @param triangles Triangles represented as sets of vertex indices.
 * @return A deduplicated vector of edges.
 */
vector<Edge> triangles_to_edges(const std::set<std::set<int>>& triangles);

inline vector<Edge> triangles_to_edges(int /*delta*/, const vector<Triangle>& triangles) {
  return triangles_to_edges(triangles);
}

inline vector<Edge> triangles_to_edges(int /*delta*/, const std::set<std::set<int>>& triangles) {
  return triangles_to_edges(triangles);
}

}

namespace Utils {

// Bring Isotopy types into Utils for convenience
using Isotopy::vector;
using Isotopy::string;
using Isotopy::pair;

// Backwards-compatible helper (master API)
std::map<std::pair<int,int>, int> get_pt2int(int delta);

// Legacy API: returns set<set<int>>
std::pair<vector<bool>, std::set<std::set<int>>> pcom_to_signs_and_triangles(const string& pcom_string);

// New API: returns vector<Triangle>
pair<vector<bool>, vector<Isotopy::Triangle>> pcom_to_signs_and_triangles_vec(const string& pcom_string);

// Backwards-compatible overload: returns set<set<int>>
std::pair<vector<bool>, std::set<std::set<int>>> pcom_to_signs_and_triangles_set(const string& pcom_string);

// New API: takes vector<Triangle>
string signs_and_triangles_to_pcom(const vector<bool>& sign_vector, const vector<Isotopy::Triangle>& triangles);
string signs_and_triangles_to_pcom(const vector<bool>& sign_vector, const vector<Isotopy::Triangle>& triangles, string origin_tag);

// Backwards-compatible overloads: take set<set<int>>
string signs_and_triangles_to_pcom(const vector<bool>& sign_vector, const std::set<std::set<int>>& triangles);
string signs_and_triangles_to_pcom(const vector<bool>& sign_vector, const std::set<std::set<int>>& triangles, string origin_tag);

}
