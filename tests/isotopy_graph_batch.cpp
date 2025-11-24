#include "catch.hpp"
#include "node.hpp"
#include "isotopy_graph.h"
#include <vector>
#include <set>
#include <array>

#include <fstream>
#include <sstream>

#include <regex>
#include <cstdio>
#include <string>
#include <string_view>
#include <iostream>



struct TestCaseData {
    int p_regions;
    int n_regions;
    int case_number;
    std::vector<std::array<int, 3>> triangulation;
    std::set<std::set<int>> triangulation_vec;
    std::string signs;
    std::vector<bool> signs_vec;
};

TestCaseData parse_test_case(const std::string& line) {
    TestCaseData data;
    std::smatch match;

    std::regex re(R"((\d+)\s+Even:\s*(\d+);\s*Odd:\s*(\d+);\s*Triangulation:\s*(\[\[.*?\]\]);\s*Signs:\s*([01]+);\s*Tree:\s*(\[.*\]))");
    if (std::regex_search(line, match, re)) {
        data.p_regions = std::stoi(match[2]);
        data.n_regions = std::stoi(match[3]);
        data.case_number = std::stoi(match[1]);
        data.signs = match[5];
        data.signs_vec.clear();
        for (char c : data.signs) data.signs_vec.push_back(c == '1');


        // Parse triangulation
        std::string tri_str = match[4];
        std::regex tri_re(R"(\[(\d+),(\d+),(\d+)\])");
        auto tri_begin = std::sregex_iterator(tri_str.begin(), tri_str.end(), tri_re);
        auto tri_end = std::sregex_iterator();
        for (auto it = tri_begin; it != tri_end; ++it) {
            data.triangulation.push_back({std::stoi((*it)[1]), std::stoi((*it)[2]), std::stoi((*it)[3])});
        }

        data.triangulation_vec.clear();
        for (const auto& tri : data.triangulation) {
          data.triangulation_vec.insert(std::set<int>{tri[0], tri[1], tri[2]});
        }


    }
    return data;
}
TEST_CASE("Isotopy::Graph batch test from tree.txt", "[isotopy_graph]") {
  std::cout << "Reading test cases from file tree.txt\n";
  std::ifstream infile("tests/tree.txt");
  REQUIRE(infile);
  // Count total lines
  std::istreambuf_iterator<char> begin(infile), end;
  int total_lines = std::count(begin, end, '\n');
  infile.clear();
  infile.seekg(0, std::ios::beg);
  int line_count = 0;
  std::string line;

  while (std::getline(infile, line)) {
      ++line_count;
      auto data = parse_test_case(line);
      int delta = 8;
      if (line_count % 500 == 0 || line_count == total_lines) {
        std::cout << "Processed " << line_count << " / " << total_lines << " test cases from tree.txt\n";
      }
    Isotopy::Graph graph(delta, data.signs_vec, data.triangulation_vec);
    graph.isotopy_type();

    int expected_p = data.p_regions;
    int expected_n = data.n_regions;

    if (graph.p_regions != expected_p || graph.n_regions != expected_n) {
      std::cout << "Discrepancy in case " << data.case_number << ": Expected (P,N)=(" << expected_p << "," << expected_n << "), Got (P,N)=(" << graph.p_regions << "," << graph.n_regions << ")\n";
    }

    REQUIRE(graph.p_regions == expected_p);
    REQUIRE(graph.n_regions == expected_n);

    std::vector<Isotopy::Triangle> triangles(data.triangulation.begin(), data.triangulation.end());
    Isotopy::Graph triangle_graph(delta, data.signs_vec, triangles);
    triangle_graph.isotopy_type();
    REQUIRE(triangle_graph.even_regions() == expected_p);
    REQUIRE(triangle_graph.odd_regions() == expected_n);
    REQUIRE(triangle_graph.viro_notation() == graph.viro_notation());
  }
}
TEST_CASE("Isotopy::Graph batch test from YAML file", "[isotopy_graph][yaml]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);
  int line_count = 0;
  int total_lines = doc.size();

  for (auto& node : doc) {
    ++line_count;

    if (line_count % 500 == 0 || line_count == total_lines) {
      std::cout << "Processed " << line_count << " / " << total_lines << " test cases from isotopy_tests.yaml\n";
    }
    int delta = node.at("degree").get_value<int>();
    std::vector<bool> signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
    std::set<std::set<int>> triangulation_vec = node.at("triangulation").get_value<std::set<std::set<int>>>();
    int expected_p = node.at("even").get_value<int>(); 
    int expected_n = node.at("odd").get_value<int>();
    size_t nverts = Isotopy::num_vertices(delta);
    if (signs_vec.size() != nverts) {
      continue; // Skip invalid test case
    }
    std::string viro = node.at("viro").get_value<std::string>();

    if (triangulation_vec.size() != static_cast<size_t>(delta * delta)) {
      std::cout << "Skipping YAML case " << line_count << " (delta=" << delta
                << ") due to invalid triangulation size: "
                << triangulation_vec.size() << "\n";
      continue;
    }
    Isotopy::Graph graph(delta, signs_vec, triangulation_vec);
    graph.isotopy_type();
    REQUIRE(graph.even_regions() == expected_p);
    REQUIRE(graph.odd_regions() == expected_n);
    REQUIRE(graph.viro_notation() == viro);
    /*
    auto sign_triangles_pair = Utils::pcom_to_signs_and_triangles(pcom);
    REQUIRE(sign_triangles_pair.first == sign);
    REQUIRE(sign_triangles_pair.second == triangles);
    std::string pcom_created =  Utils::signs_and_triangles_to_pcom(sign_triangles_pair.first, sign_triangles_pair.second);
    auto sign_triangles_pair_created = Utils::pcom_to_signs_and_triangles(pcom_created);
    REQUIRE(sign_triangles_pair_created.first == sign);
    REQUIRE(sign_triangles_pair_created.second == triangles);
    */
    std::string pcom = Utils::signs_and_triangles_to_pcom(signs_vec, triangulation_vec);
    std::pair<std::vector<bool>, std::set<std::set<int>>> sign_triangles_pair = Utils::pcom_to_signs_and_triangles(pcom);
    REQUIRE(sign_triangles_pair.first == signs_vec);
    REQUIRE(sign_triangles_pair.second == triangulation_vec);

    std::vector<Isotopy::Triangle> tri_vec;
    tri_vec.reserve(triangulation_vec.size());
    for (const auto& tri : triangulation_vec) {
      std::vector<int> sorted_tri(tri.begin(), tri.end());
      tri_vec.push_back({sorted_tri[0], sorted_tri[1], sorted_tri[2]});
    }
    Isotopy::Graph triangle_graph(delta, signs_vec, tri_vec);
    triangle_graph.isotopy_type();
    REQUIRE(triangle_graph.even_regions() == expected_p);
    REQUIRE(triangle_graph.odd_regions() == expected_n);
    REQUIRE(triangle_graph.viro_notation() == viro);


  }
}



TEST_CASE("Isotopy::Graph batch test from mcurves.txt.xz", "[isotopy_graph]") {
  std::cout << "Reading test cases from compressed file mcurves.txt.xz\n";
  FILE* pipe = popen("xz -dc tests/mcurves.txt.xz", "r");
  REQUIRE(pipe != nullptr);

  // Count total lines
  FILE* count_pipe = popen("xz -dc tests/mcurves.txt.xz | wc -l", "r");
  int total_lines = 0;
  fscanf(count_pipe, "%d", &total_lines);
  pclose(count_pipe);

  int line_count = 0;
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), pipe)) {
    ++line_count;
    std::string line(buffer);
    auto data = parse_test_case(line);
    int delta = 8;
    if (data.case_number == 1202044) {
      continue;
    }
    if (line_count % 500 == 0 || line_count == total_lines) {
      std::cout << "Processed " << line_count << " / " << total_lines << " test cases from mcurves.txt.xz\n";
    }
    Isotopy::Graph graph(delta, data.signs_vec, data.triangulation_vec);
    graph.isotopy_type();

    int expected_p = data.p_regions;
    int expected_n = data.n_regions;

    REQUIRE(graph.p_regions == expected_p);
    REQUIRE(graph.n_regions == expected_n);

    std::vector<Isotopy::Triangle> triangles(data.triangulation.begin(), data.triangulation.end());
    Isotopy::Graph triangle_graph(delta, data.signs_vec, triangles);
    triangle_graph.isotopy_type();
    REQUIRE(triangle_graph.even_regions() == expected_p);
    REQUIRE(triangle_graph.odd_regions() == expected_n);
    REQUIRE(triangle_graph.viro_notation() == graph.viro_notation());
  }
  pclose(pipe);
}
