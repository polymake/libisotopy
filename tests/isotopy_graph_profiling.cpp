#include "isotopy_graph.h"
#include <vector>
#include <string>
#include <cstdlib>

// Usage: ./isotopy_graph_profiling [mode] [iterations]
// Modes: triangles (default), edges, phases-triangles, phases-edges
// Iterations: default 1000

int main(int argc, char* argv[]) {
  std::string mode = (argc > 1) ? argv[1] : "triangles";
  int iterations = (argc > 2) ? std::atoi(argv[2]) : 1000;
  int delta = 8;

  // Test data (delta=8 case)
  std::vector<bool> sign {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1};

  std::vector<Isotopy::Triangle> triangles {
    {0, 1, 9}, {1, 2, 11}, {1, 9, 10}, {1, 10, 19}, {1, 11, 19}, {2, 3, 11},
    {3, 4, 12}, {3, 11, 19}, {3, 12, 19}, {4, 5, 12}, {5, 6, 13}, {5, 12, 19},
    {5, 13, 19}, {6, 7, 13}, {7, 8, 14}, {7, 13, 19}, {7, 14, 19}, {8, 14, 15},
    {8, 15, 16}, {9, 10, 19}, {9, 17, 18}, {9, 18, 19}, {14, 15, 19}, {15, 16, 19},
    {16, 19, 20}, {16, 20, 21}, {16, 21, 22}, {16, 22, 23}, {17, 18, 24}, {18, 19, 24},
    {19, 20, 29}, {19, 24, 25}, {19, 25, 35}, {19, 26, 38}, {19, 26, 43}, {19, 27, 28},
    {19, 27, 33}, {19, 28, 29}, {19, 31, 35}, {19, 31, 42}, {19, 33, 38}, {19, 36, 40},
    {19, 36, 42}, {19, 40, 43}, {20, 21, 29}, {21, 22, 29}, {22, 23, 29}, {24, 25, 30},
    {25, 30, 35}, {26, 32, 38}, {26, 32, 43}, {27, 28, 33}, {28, 29, 34}, {28, 33, 34},
    {31, 35, 39}, {31, 39, 42}, {32, 37, 38}, {32, 37, 43}, {33, 34, 38}, {36, 40, 44},
    {36, 42, 44}, {37, 38, 41}, {37, 41, 43}, {40, 43, 44}
  };

  std::vector<Isotopy::Edge> edges = Isotopy::triangles_to_edges(triangles);

  for (int i = 0; i < iterations; ++i) {
    if (mode == "triangles") {
      Isotopy::Graph graph(delta, sign, triangles);
      graph.isotopy_type();
      volatile auto notation = graph.viro_notation();

    } else if (mode == "edges") {
      Isotopy::Graph graph(delta, sign, edges);
      graph.isotopy_type();
      volatile auto notation = graph.viro_notation();

    } else if (mode == "phases-triangles") {
      Isotopy::Graph graph(delta);
      graph.initialize(sign);
      graph.process_triangles(triangles);
      graph.isotopy_type();
      volatile auto notation = graph.viro_notation();

    } else if (mode == "phases-edges") {
      Isotopy::Graph graph(delta);
      graph.initialize(sign);
      graph.process_edges(edges);
      graph.isotopy_type();
      volatile auto notation = graph.viro_notation();
    }
  }

  return 0;
}
