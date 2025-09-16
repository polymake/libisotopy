#include "catch.hpp"
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

    REQUIRE(graph.p_regions == expected_p);
    REQUIRE(graph.n_regions == expected_n);
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
  }
  pclose(pipe);
}

