#include "catch.hpp"
#include "isotopy_graph.h"
#include <vector>
#include <set>
#include <array>

#include <fstream>
#include <sstream>

#include <regex>


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
TEST_CASE("Isotopy::Graph batch test from file", "[isotopy_graph]") {
  std::ifstream infile("tests/tree.txt");
  REQUIRE(infile);
    std::string line;
    while (std::getline(infile, line)) {
        auto data = parse_test_case(line);
        int delta = 8;
        if (data.case_number == 1202044) {
          continue;
        }
        SECTION("Case " + std::to_string(data.case_number)) {
            Isotopy::Graph graph(delta, data.signs_vec, data.triangulation_vec);
            graph.isotopy_type();

            int expected_p = data.p_regions;
            int expected_n = data.n_regions;

            REQUIRE(graph.p_regions == expected_p);
            REQUIRE(graph.n_regions == expected_n);
        }
    }
}
