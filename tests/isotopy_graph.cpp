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
        for (size_t i = 0; i < graph.component_adjacency.size(); ++i) {
            std::cout << "Component " << i << ": ";
            for (size_t ele = 0; ele < graph.component.size(); ++ele) {
                if (graph.component[ele] == static_cast<int>(i)) {
                    std::cout << ele << " ";
                }
            }
            std::cout << "\n";
        }
    }

    if (!graph.nbs.empty()) {
        std::cout << "\nNBS:\n";
        for (size_t i = 0; i < graph.nbs.size(); ++i) {
            std::cout << "Side point "<< i<< " is point: " << graph.side_points[i] << ": ";
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
    if (graph.root_region != -1) {
        std::cout << "Root Region: " << graph.root_region << "\n";
    }
    if (graph.region.size() > 0 && graph.component.size() > 0) {
        std::cout << "\nRegions:\n";
        for (size_t i = 0; i < graph.region.size(); ++i) {
            std::cout << "Region " << i << ": ";
            for (size_t vert = 0; vert < graph.component.size(); ++vert) {
                int comp = graph.component[vert];
                if (graph.region[comp] == static_cast<int>(i)) {
                    std::cout << vert << " ";
                }
            }
            std::cout << "\n";
        }
    }
    if (graph.region_adjacency.size() > 0) {
        std::cout << "\nRegion Adjacency:\n";
        for (size_t i = 0; i < graph.region_adjacency.size(); ++i) {
            std::cout << "Region " << i << ": ";
            for (int adj : graph.region_adjacency[i]) {
                std::cout << adj << " ";
            }
            std::cout << "\n";
        }
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
        graph.isotopy_root();
        REQUIRE(graph.root == 1);
        graph.isotopy_type();
        REQUIRE(graph.region_sign == std::vector<bool>({true, false}));
        REQUIRE(graph.p_regions == 1);
        REQUIRE(graph.n_regions == 0);
    }
    SECTION("Viro notation") {
        auto viro = graph.viro_notation();
        REQUIRE(viro == "<1>");
    }
    SECTION("shorthand") {
        auto shorthand = Utils::shorthand(delta, sign, edges);
        REQUIRE(shorthand == "Ahqw~");
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
    std::string expected_viro = "<4v1<1>>";
    REQUIRE(graph.viro_notation() == expected_viro);
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
    REQUIRE(graph.viro_notation() == "<3v1<1>>");

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

    REQUIRE(graph.viro_notation() == "<2v1<2>>");
}

TEST_CASE("Isotopy::Graph Test Case 29", "[isotopy_graph]") {
    int delta = 8;
    std::vector<bool> sign {1,1,0,0,1,0,0,1,1,0,1,0,1,0,0,1,1,0,0,0,0,1,1,0,1,1,1,0,0,1,0,0,1,1,1,0,1,0,1,1,0,1,1,1,0};

    std::set<std::set<int>> triangles {{0,9,10},{0,1,10},{9,10,17},{1,2,10},{10,17,18},{2,10,11},{10,18,19},{10,11,19},{17,18,25},{2,11,12},{18,19,25},{11,12,19},{17,24,25},{2,3,12},{24,25,30},{19,25,26},{12,19,20},{3,4,12},{19,26,27},{19,20,27},{25,26,32},{12,20,21},{25,30,31},{4,12,13},{25,31,32},{26,27,32},{20,21,27},{12,13,21},{30,31,36},{4,13,14},{31,32,36},{13,14,21},{27,32,33},{21,27,28},{30,35,36},{4,5,14},{27,33,34},{27,28,34},{32,36,37},{32,33,38},{21,28,29},{14,21,22},{35,36,39},{32,37,38},{21,22,29},{5,6,14},{33,34,38},{28,29,34},{36,37,41},{14,22,23},{37,38,41},{22,23,29},{36,39,40},{6,14,15},{36,40,41},{14,15,23},{39,40,43},{6,15,16},{40,41,43},{15,16,23},{39,42,43},{6,7,16},{42,43,44},{7,8,16}};
    Isotopy::Graph graph(delta, sign, triangles);


    int expected_p = 5;
    int expected_n = 3;

    REQUIRE(graph.even_regions() == expected_p);
    REQUIRE(graph.odd_regions() == expected_n);
    REQUIRE(graph.viro_notation() == "<3v1<1>v1<2>>");

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
    REQUIRE(graph.viro_notation() == "<8v1<6v1<3>>>");

}

TEST_CASE("Isotopy::Graph Sebastian's example", "[isotopy_graph]") {

    int delta = 6;
    std::vector<bool> sign {0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,0,0,1,0,1,0,1,1,1,1,0,1,1};
    std::set<std::set<int>> triangles {{5,6,12},{5,11,12},{4,5,11},{11,12,17},{4,10,11},{11,16,17},{3,4,10},{10,11,16},{16,17,21},{2,3,10},{2,9,10},{9,10,16},{1,2,9},{9,16,21},{1,9,21},{1,15,21},{1,8,15},{1,7,8},{7,8,15},{0,1,7},{7,15,21},{7,14,21},{14,20,21},{14,19,20},{7,13,14},{13,14,19},{20,21,24},{20,23,24},{19,20,23},{23,24,26},{13,18,19},{19,22,23},{18,19,22},{23,25,26},{22,23,25},{25,26,27}};
    Isotopy::Graph graph(delta, sign, triangles);
    graph.isotopy_type();
    REQUIRE(graph.p_regions == 1);
    REQUIRE(graph.n_regions == 9);
    REQUIRE(graph.viro_notation() == "<1<9>>");
}
TEST_CASE("Isotopy::Graph 14679", "[isotopy_graph]") {

    int delta = 8;
    std::vector<bool> sign {1,1,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,0,0,1,1,1,1,0,1,0,0,1,0,1};
    std::set<std::set<int>> triangles {{8,15,16},{7,8,15},{15,16,23},{6,7,15},{6,14,28},{6,15,22},{15,22,34},{5,21,27},{15,23,29},{5,21,28},{5,6,14},{22,28,34},{5,14,28},{6,22,28},{3,4,13},{15,29,34},{27,28,31},{19,20,27},{4,5,13},{5,13,27},{28,34,38},{12,13,17},{28,32,33},{2,3,13},{28,33,38},{13,20,27},{13,19,20},{26,27,31},{18,19,27},{21,27,28},{2,11,12},{32,33,37},{13,18,19},{33,37,38},{2,12,13},{25,26,31},{28,31,32},{37,38,40},{11,12,17},{32,36,42},{17,26,27},{17,18,27},{31,32,35},{2,10,11},{38,41,42},{10,11,17},{32,37,40},{17,24,31},{17,25,31},{32,35,36},{1,2,9},{38,40,42},{13,17,18},{24,30,31},{30,31,35},{36,39,42},{17,25,26},{35,36,39},{0,1,9},{41,42,43},{9,10,17},{32,40,42},{2,9,10},{42,43,44}};

    Isotopy::Graph graph(delta, sign, triangles);
    graph.isotopy_type();
    REQUIRE(graph.even_regions() == 2);
    REQUIRE(graph.odd_regions() == 15);
    REQUIRE(graph.viro_notation() == "<1v1<15>>");
    //std::cout << "Shorthand: " << Utils::web_shorthand(delta, sign, triangles) << std::endl;
}


