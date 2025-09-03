#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "isotopy_graph.h"
#include <vector>
#include <set>

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
    }
}
