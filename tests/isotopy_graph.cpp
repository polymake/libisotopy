#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "isotopy_graph.h"
#include <vector>
#include <set>
#include <array>

#include <fstream>
#include <sstream>

#include <regex>

void graph_debug(const Isotopy::Graph& graph) {
    std::cout << "Adjacency List:\n";
    for (size_t i = 0; i < graph.adjacency.size(); ++i) {
        std::cout << i << ": ";
        for (int neighbor : graph.adjacency[i]) {
            std::cout << neighbor << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nSides:\n";
    for (size_t i = 0; i < graph.sides.size(); ++i) {
        std::cout << "Side point " << graph.side_points[i] << ": ";
        for (bool side : graph.sides[i]) {
            std::cout << side << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nSign Vector:\n";
    for (size_t i = 0; i < graph.sign.size(); ++i) {
        std::cout << graph.sign[i] << " ";
    }
    std::cout << "\n";

    if (!graph.component.empty()) {
        std::cout << "\nComponents:\n";
        for (size_t i = 0; i < graph.component.size(); ++i) {
            std::cout << "Vertex " << i << ": Component " << graph.component[i] << "\n";
        }
    }

    if (!graph.nbs.empty()) {
        std::cout << "\nNBS:\n";
        for (size_t i = 0; i < graph.nbs.size(); ++i) {
            std::cout << "Side point " << graph.side_points[i] << ": ";
            for (int nb : graph.nbs[i]) {
                std::cout << nb << " ";
            }
            std::cout << "\n";
        }
    }

    if (graph.root != -1) {
        std::cout << "\nIsotopy Root: " << graph.root << "\n";
    }

    if (!graph.region_sign.empty()) {
        std::cout << "\nRegion Signs:\n";
        for (size_t i = 0; i < graph.region_sign.size(); ++i) {
            std::cout << "Region " << i << ": Sign " << graph.region_sign[i] << "\n";
        }
    }

    if (graph.p_regions != -1 && graph.n_regions != -1) {
        std::cout << "\nPositive Regions: " << graph.p_regions << "\n";
        std::cout << "Negative Regions: " << graph.n_regions << "\n";
    }
}

TEST_CASE("Isotopy::Graph basic properties", "[isotopy_graph]") {
    int delta = 2;
    std::vector<bool> sign = {true, true, true, true, true, true};
    std::set<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {0, 3}, {3, 4}, {2, 4}, {1,3}, {1,4}, {4,5}, {3,5}};
    Isotopy::Graph graph(delta, sign, edges);

    size_t nverts_triang = (delta + 1) * (delta + 2) / 2;
    size_t nvert_formula = 4 * nverts_triang - 4 * (delta + 1) + 1;

    SECTION("Adjacency size") {
        REQUIRE(graph.adjacency.size() == nvert_formula);
    }

    SECTION("Adjacency details") {
        REQUIRE(graph.adjacency[0].size() == 4);
        REQUIRE(graph.adjacency[1].size() == 6);
        REQUIRE(graph.adjacency[10].size() == 4);
        REQUIRE(graph.adjacency[12].size() == 3);
        REQUIRE(graph.adjacency[9].size() == 6);
    }

    SECTION("Side points and sides") {
        REQUIRE(graph.sides.size() == 8);
        REQUIRE(graph.sides[0] == std::vector<bool>({true, false, false, true}));
        REQUIRE(graph.sides[1] == std::vector<bool>({true, false, false, false}));
        REQUIRE(graph.sides[2] == std::vector<bool>({true, true, false, false}));
        REQUIRE(graph.sides[3] == std::vector<bool>({false, true, false, false}));
        REQUIRE(graph.sides[4] == std::vector<bool>({false, true, true, false}));
        REQUIRE(graph.sides[5] == std::vector<bool>({false, false, true, false}));
        REQUIRE(graph.sides[6] == std::vector<bool>({false, false, true, true}));
        REQUIRE(graph.sides[7] == std::vector<bool>({false, false, false, true}));
        REQUIRE(graph.side_points == std::vector<int>({2, 4, 5, 8, 7, 10, 12, 11}));
        REQUIRE(graph.side_points.size() == 8);
    }

    SECTION("Sign vector") {
        REQUIRE(graph.sign == std::vector<bool>({true, true, true, true, true, true, false, true, false, false, true, false, true}));
    }

    SECTION("Connected components") {
        graph.connected_components();
        REQUIRE(graph.component.size() == graph.adjacency.size());
        REQUIRE(graph.component_adjacency.size() == 3);
        REQUIRE(graph.component == std::vector<int>({0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 2, 1, 2}));
    }

    SECTION("Pre isotopy root and nbs") {
        graph.connected_components();
        graph.pre_isotopy_root();
        REQUIRE(graph.nbs.size() == graph.side_points.size());
        REQUIRE(graph.nbs[0] == std::set<int>({1, 2}));
        REQUIRE(graph.nbs[1] == std::set<int>({0, 2}));
        REQUIRE(graph.nbs[2] == std::set<int>({0, 1}));
        REQUIRE(graph.nbs[3] == std::set<int>({7}));
        REQUIRE(graph.nbs[4] == std::set<int>({5, 6}));
        REQUIRE(graph.nbs[5] == std::set<int>({4, 6}));
        REQUIRE(graph.nbs[6] == std::set<int>({4, 5}));
        REQUIRE(graph.nbs[7] == std::set<int>({3}));
    }

    SECTION("Isotopy root") {
        graph.connected_components();
        graph.pre_isotopy_root();
        graph.isotopy_root();
        REQUIRE(graph.root == 1);
        graph.calculate_regions();
        graph.isotopy_type();
        REQUIRE(graph.region_sign == std::vector<bool>({true, false}));
        REQUIRE(graph.p_regions == 1);
        REQUIRE(graph.n_regions == 0);
    }
}

TEST_CASE("Isotopy::Graph from triangulation ", "[isotopy_graph]") {
    int delta = 8;
    std::vector<bool> sign {0,1,0,0,0,1,0,0,1,1,1,0,0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,0,0,1};
    std::set<std::set<int>> triangles {{0,9,10},{0,1,10},{9,10,17},{1,2,10},{10,17,18},{2,10,11},{10,18,19},{10,11,19},{17,18,25},{2,11,12},{18,19,25},{11,12,19},{17,24,25},{2,3,12},{24,25,30},{19,25,26},{12,19,20},{3,4,12},{19,26,27},{19,20,27},{25,26,32},{12,20,21},{25,30,31},{4,12,13},{25,31,32},{26,27,32},{20,21,27},{12,13,21},{30,31,36},{4,13,14},{31,32,36},{13,14,21},{27,32,33},{21,27,28},{30,35,36},{4,5,14},{27,33,34},{27,28,34},{32,36,37},{32,33,38},{21,28,29},{14,21,22},{35,36,39},{32,37,38},{21,22,29},{5,6,14},{33,34,38},{28,29,34},{36,37,41},{14,22,23},{37,38,41},{22,23,29},{36,39,40},{6,14,15},{36,40,41},{14,15,23},{39,40,43},{6,15,16},{40,41,43},{15,16,23},{39,42,43},{6,7,16},{42,43,44},{7,8,16}};

    Isotopy::Graph graph(delta, sign, triangles);
    graph.isotopy_type();
    REQUIRE(graph.p_regions == 5);
    REQUIRE(graph.n_regions == 1);
}

TEST_CASE("Isotopy::Graph Test Case 9", "[isotopy_graph]") {
    int delta = 8;
    std::vector<bool> sign {0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,0,0,1,0,1,0,1,0,0,1,0,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,0,1,1,1};
    std::set<std::set<int>> triangles {{0,9,10},{0,1,10},{9,10,17},{1,2,10},{10,17,18},{2,10,11},{10,18,19},{10,11,19},{17,18,25},{2,11,12},{18,19,25},{11,12,19},{17,24,25},{2,3,12},{24,25,30},{19,25,26},{12,19,20},{3,4,12},{19,26,27},{19,20,27},{25,26,32},{12,20,21},{25,30,31},{4,12,13},{25,31,32},{26,27,32},{20,21,27},{12,13,21},{30,31,36},{4,13,14},{31,32,36},{13,14,21},{27,32,33},{21,27,28},{30,35,36},{4,5,14},{27,33,34},{27,28,34},{32,36,37},{32,33,38},{21,28,29},{14,21,22},{35,36,39},{32,37,38},{21,22,29},{5,6,14},{33,34,38},{28,29,34},{36,37,41},{14,22,23},{37,38,41},{22,23,29},{36,39,40},{6,14,15},{36,40,41},{14,15,23},{39,40,43},{6,15,16},{40,41,43},{15,16,23},{39,42,43},{6,7,16},{42,43,44},{7,8,16}};
    Isotopy::Graph graph(delta, sign, triangles);
    int expected_p = 4;
    int expected_n = 1;
    graph.isotopy_type();


    REQUIRE(graph.p_regions == expected_p);
    REQUIRE(graph.n_regions == expected_n);


}

TEST_CASE("Isotopy::Graph Test Case 34", "[isotopy_graph]") {
    int delta = 8;
    std::vector<bool> sign {1,1,1,0,0,1,0,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,1,0,0,1,0,0,1,0,1,0,1,0,1,1,0,0,1,1,0,1,1,0};
    std::set<std::set<int>> triangles {{0,9,10},{0,1,10},{9,10,17},{1,2,10},{10,17,18},{2,10,11},{10,18,19},{10,11,19},{17,18,25},{2,11,12},{18,19,25},{11,12,19},{17,24,25},{2,3,12},{24,25,30},{19,25,26},{12,19,20},{3,4,12},{19,26,27},{19,20,27},{25,26,32},{12,20,21},{25,30,31},{4,12,13},{25,31,32},{26,27,32},{20,21,27},{12,13,21},{30,31,36},{4,13,14},{31,32,36},{13,14,21},{27,32,33},{21,27,28},{30,35,36},{4,5,14},{27,33,34},{27,28,34},{32,36,37},{32,33,38},{21,28,29},{14,21,22},{35,36,39},{32,37,38},{21,22,29},{5,6,14},{33,34,38},{28,29,34},{36,37,41},{14,22,23},{37,38,41},{22,23,29},{36,39,40},{6,14,15},{36,40,41},{14,15,23},{39,40,43},{6,15,16},{40,41,43},{15,16,23},{39,42,43},{6,7,16},{42,43,44},{7,8,16}};

    Isotopy::Graph graph = Isotopy::Graph(delta, sign, triangles);


    int expected_p = 3;
    int expected_n = 2;
    graph.isotopy_type();

    REQUIRE(graph.p_regions == expected_p);
    REQUIRE(graph.n_regions == expected_n);
}

TEST_CASE("Isotopy::Graph Test Case 29", "[isotopy_graph]") {
    int delta = 8;
    std::vector<bool> sign {1,1,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,1,0,1,1,1,0,0,1,0,0,1,1,1,0,1,0,1,1,0,1,1,1,0};

    std::set<std::set<int>> triangles {{0,9,10},{0,1,10},{9,10,17},{1,2,10},{10,17,18},{2,10,11},{10,18,19},{10,11,19},{17,18,25},{2,11,12},{18,19,25},{11,12,19},{17,24,25},{2,3,12},{24,25,30},{19,25,26},{12,19,20},{3,4,12},{19,26,27},{19,20,27},{25,26,32},{12,20,21},{25,30,31},{4,12,13},{25,31,32},{26,27,32},{20,21,27},{12,13,21},{30,31,36},{4,13,14},{31,32,36},{13,14,21},{27,32,33},{21,27,28},{30,35,36},{4,5,14},{27,33,34},{27,28,34},{32,36,37},{32,33,38},{21,28,29},{14,21,22},{35,36,39},{32,37,38},{21,22,29},{5,6,14},{33,34,38},{28,29,34},{36,37,41},{14,22,23},{37,38,41},{22,23,29},{36,39,40},{6,14,15},{36,40,41},{14,15,23},{39,40,43},{6,15,16},{40,41,43},{15,16,23},{39,42,43},{6,7,16},{42,43,44},{7,8,16}};
    Isotopy::Graph graph(delta, sign, triangles);


    int expected_p = 5;
    int expected_n = 3;

    graph.isotopy_type();

    REQUIRE(graph.p_regions == expected_p);
    REQUIRE(graph.n_regions == expected_n);

}

TEST_CASE("Isotopy::Graph Test Case 1202044", "[isotopy_graph]") {
    int delta = 8;
    std::vector<bool> sign {0,0,1,1,1,0,0,1,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1,0,0,0,1,1,1,1,0,1,1,1,1,0,0,0};
    std::set<std::set<int>> triangles {{0,1,9},{1,2,9},{2,3,9},{3,4,11},{3,9,10},{3,10,11},{4,5,12},{4,11,18},{4,12,24},{4,18,24},{5,6,13},{5,12,19},{5,13,30},{5,19,25},{5,25,30},{6,7,14},{6,13,20},{6,14,35},{6,20,26},{6,26,31},{6,31,35},{7,8,15},{7,14,21},{7,15,39},{7,21,27},{7,27,32},{7,32,36},{7,36,39},{8,15,22},{8,16,42},{8,22,28},{8,28,33},{8,33,37},{8,37,40},{8,40,42},{9,10,17},{10,11,17},{11,17,18},{12,19,24},{13,20,30},{14,21,35},{15,22,39},{16,23,42},{17,18,24},{19,24,25},{20,26,30},{21,27,35},{22,28,39},{23,29,42},{24,25,30},{26,30,31},{27,32,35},{28,33,39},{29,34,42},{30,31,35},{32,35,36},{33,37,39},{34,38,42},{35,36,39},{37,39,40},{38,41,42},{39,40,42},{41,42,43},{42,43,44}};
    Isotopy::Graph graph(delta, sign, triangles);

    int expected_p = 12; 
    int expected_n = 7;

    graph.isotopy_type();

    //graph_debug(graph);
    REQUIRE(graph.p_regions == expected_p);
    REQUIRE(graph.n_regions == expected_n);

}


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
