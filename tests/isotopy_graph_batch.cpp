#include "catch.hpp"
#include "node.hpp"
#include "isotopy_graph.h"
#include <vector>
#include <set>
#include <array>
#include <fstream>
#include <string>
#include <iostream>

TEST_CASE("Batch Tests: isotopy_tests.yaml", "[isotopy_graph][batch]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  if (!ifs.is_open()) {
    WARN("Could not open tests/isotopy_tests.yaml. Skipping batch tests.");
    return;
  }
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);
  int line_count = 0;
  int total_lines = doc.size();

  for (auto& node : doc) {
    ++line_count;

    if (line_count % 500 == 0 || line_count == total_lines) {
      std::cout << "Processed " << line_count << " / " << total_lines << " test cases from isotopy_tests.yaml\n";
    }

    // Parse test case with clear error reporting
    std::string case_id = "case #" + std::to_string(line_count);
    try {
      if (node.contains("comment")) {
        case_id += " (" + node.at("comment").get_value<std::string>() + ")";
      }
    } catch (...) {}

    int delta = 0;
    try {
      delta = node.at("degree").get_value<int>();
    } catch (const std::exception& e) {
      FAIL("Missing or invalid 'degree' field in " << case_id << ": " << e.what());
    }

    std::vector<bool> signs_vec;
    try {
      signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
    } catch (const std::exception& e) {
      FAIL("Missing or invalid 'polarisation' field in " << case_id << ": " << e.what());
    }

    std::set<std::set<int>> triangulation_set;
    try {
      triangulation_set = node.at("triangulation").get_value<std::set<std::set<int>>>();
    } catch (const std::exception& e) {
      FAIL("Missing or invalid 'triangulation' field in " << case_id << ": " << e.what());
    }

    // Parse edges if available, otherwise derive
    std::vector<std::pair<int,int>> edges_vec;
    bool edges_from_yaml = false;
    if (node.contains("edges")) {
        try {
            edges_vec = node.at("edges").get_value<std::vector<std::pair<int,int>>>();
            edges_from_yaml = true;
        } catch (const std::exception& e) {
            WARN("Failed to parse 'edges' in " << case_id << ": " << e.what() << ". Deriving from triangulation.");
        }
    }

    int expected_p, expected_n;
    std::string expected_viro;
    try {
      expected_p = node.at("even").get_value<int>();
      expected_n = node.at("odd").get_value<int>();
      if (node.contains("viro")) {
        expected_viro = node.at("viro").get_value<std::string>();
      } else {
        expected_viro = "UNKNOWN";
      }
    } catch (const std::exception& e) {
      FAIL("Missing or invalid even/odd fields in " << case_id << ": " << e.what());
    }

    // Validate data before constructing graph
    size_t nverts = Isotopy::num_vertices(delta);
    size_t ntriangles = delta * delta;

    if (signs_vec.size() != nverts) {
      WARN("Skipping " << case_id << ": wrong number of signs (" << signs_vec.size() << " != " << nverts << ")");
      continue;
    }

    if (triangulation_set.size() != ntriangles) {
      WARN("Skipping " << case_id << ": wrong number of triangles (" << triangulation_set.size() << " != " << ntriangles << ")");
      continue;
    }

    // --- Constructor 1: set<set<int>> (Triangulation) ---
    {
        Isotopy::Graph graph(delta, signs_vec, triangulation_set);
        graph.isotopy_type();
        CHECK(graph.even_regions() == expected_p);
        CHECK(graph.odd_regions() == expected_n);
        if (expected_viro != "UNKNOWN") {
            CHECK(graph.viro_notation() == expected_viro);
        }
    }

    // --- Constructor 2: vector<Triangle> ---
    {
        std::vector<Isotopy::Triangle> triangles_vec;
        triangles_vec.reserve(triangulation_set.size());
        for(const auto& t : triangulation_set) {
            std::vector<int> tv(t.begin(), t.end());
            triangles_vec.push_back({tv[0], tv[1], tv[2]});
        }
        Isotopy::Graph graph(delta, signs_vec, triangles_vec);
        graph.isotopy_type();
        CHECK(graph.even_regions() == expected_p);
        CHECK(graph.odd_regions() == expected_n);
        if (expected_viro != "UNKNOWN") {
            CHECK(graph.viro_notation() == expected_viro);
        }

        // If edges weren't in YAML, derive them now for next tests
        if (!edges_from_yaml) {
             edges_vec = Isotopy::triangles_to_edges(triangles_vec);
        }
    }

    // --- Constructor 3: vector<Edge> ---
    {
        Isotopy::Graph graph(delta, signs_vec, edges_vec);
        graph.isotopy_type();
        CHECK(graph.even_regions() == expected_p);
        CHECK(graph.odd_regions() == expected_n);
        if (expected_viro != "UNKNOWN") {
            CHECK(graph.viro_notation() == expected_viro);
        }
    }

    // --- Constructor 4: set<Edge> ---
    {
        std::set<std::pair<int,int>> edges_set(edges_vec.begin(), edges_vec.end());
        Isotopy::Graph graph(delta, signs_vec, edges_set);
        graph.isotopy_type();
        CHECK(graph.even_regions() == expected_p);
        CHECK(graph.odd_regions() == expected_n);
        if (expected_viro != "UNKNOWN") {
            CHECK(graph.viro_notation() == expected_viro);
        }
    }
  }
}
