#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "isotopy_graph.h"
#include <vector>
#include <set>
#include <array>

#include <fstream>
#include <sstream>

#include <regex>


TEST_CASE("Isotopy::Graph from triangulation  (Harnack 8)", "[isotopy_graph]") {
  int delta = 8;
  std::vector<bool> sign {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  std::set<std::set<int>> triangles = {{1,0,9},{1,9,10},{2,1,10},{2,10,11},{3,2,11},{3,11,12},{4,3,12},{4,12,13},{5,4,13},{5,13,14},{6,5,14},{6,14,15},{7,6,15},{7,15,16},{8,7,16},{10,9,17},{10,17,18},{11,10,18},{11,18,19},{12,11,19},{12,19,20},{13,12,20},{13,20,21},{14,13,21},{14,21,22},{15,14,22},{15,22,23},{16,15,23},{18,17,24},{18,24,25},{19,18,25},{19,25,26},{20,19,26},{20,26,27},{21,20,27},{21,27,28},{22,21,28},{22,28,29},{23,22,29},{25,24,30},{25,30,31},{26,25,31},{26,31,32},{27,26,32},{27,32,33},{28,27,33},{28,33,34},{29,28,34},{31,30,35},{31,35,36},{32,31,36},{32,36,37},{33,32,37},{33,37,38},{34,33,38},{36,35,39},{36,39,40},{37,36,40},{37,40,41},{38,37,41},{40,39,42},{40,42,43},{41,40,43},{43,42,44}};

  Isotopy::Graph graph(delta, sign, triangles);
  graph.isotopy_type();
  std::string expected_viro = "<1<1<1<1>>>>";
  REQUIRE(graph.viro_notation() == expected_viro);
  std::map<std::pair<int,int>, int> pt_map = Utils::get_pt2int(delta);
  REQUIRE(graph.component[pt_map.at(std::make_pair(0,0))] == graph.component[pt_map.at(std::make_pair(0,1))]);
  REQUIRE(graph.component[pt_map.at(std::make_pair(0,0))] == graph.component[pt_map.at(std::make_pair(1,0))]);
  REQUIRE(graph.component[pt_map.at(std::make_pair(0,-1))] == graph.component[pt_map.at(std::make_pair(1,-1))]);
  //std::string pcom = Utils::signs_and_triangles_to_pcom(sign,triangles);
  //REQUIRE(Utils::pcom_to_signs_and_triangles(pcom) == std::make_pair(sign,triangles)); //Maybe sort triangles

}

TEST_CASE("Isotopy::Graph from triangulation  (Loading from pcom Test)", "[isotopy_graph]") {
    std::string pcom = "{\"TYPE\":        \"<13v1<6>v1<1>>\",\"PATCHWORK\":[{\"SIGNS\":\n[true,true,true,true,true,true,true,true,false,true,false,true,true,true,true,true,\n true,true,true,true,true,true,true,false,true,false,true,true,true,true,true,true,true,true,false,true,true,true,true,true,false,false,true,true]}],\"_ns\":{\"polymake\":[\"https://polymake.org\",\"4.13\"]},\"DUAL_SUBDIVISION\":{\"WEIGHTS\":[\"113\",\"321\",\"530\",\"740\",\"951\",\"1163\",\"1376\",\"1590\",\"1805\",\"96\",\"35\",\"154\",\"355\",\"559\",\"764\",\"970\",\"1177\",\"80\",\"26\",\"9\",\"0\",\"196\",\"395\",\"598\",\"65\",\"18\",\"4\",\"4\",\"18\",\"40\",\"51\",\"11\",\"0\",\"9\",\"28\",\"38\",\"5\",\"1\",\"17\",\"26\",\"2\",\"7\",\"15\",\"0\",\"5\"],\"MAXIMAL_CELLS\":[[0,9,10],[10,18,19],[20,27,28],[0,12,13],[0,13,14],[0,14,15],[0,15,16],[0,1,16],[2,3,16],[1,2,16],[7,8,16],[6,7,16],[5,6,16],[4,5,16],[3,4,16],[9,10,17],[13,14,23],[12,13,23],\n[0,\n              10,\n        20],[10,19,20],[28,29,34],[0,11,22],[11,21,22],[0,12,22],[12,22,23],[15,16,23],[14,15,23],[10,17,24],[0,20,29],[0,11,29],[20,28,29],[18,19,25],[10,24,30],[19,20,26],[22,23,29],[11,21,29],[21,22,29],[20,27,37],[27,28,33],[28,34,38],[10,30,35],[19,25,31],[28,33,41],[20,26,32],[19,26,36],[20,37,43],[27,33,37],[20,32,43],[26,32,36],[19,31,36],[28,38,41],[37,41,43],[33,37,41],[10,35,39],[32,36,40],[32,40,43],[10,39,42],[10,42,44],[40,43,44],[36,40,44],[31,36,44],[25,31,44],[18,25,44],[10,18,44],{\"cols\":45}]},\"PURE\":true,\"_type\":\"tropical::Hypersurface<Min>\",\"_id\":\"15-7_nn13-c1-6-c1-1\",\"_attrs\":{\"TYPE\":{\"attachment\":true}}}";
    //std::string pcom("SIGNS");
   std::vector<bool> sign {1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,0,0,1,1};
   std::set<std::set<int>> triangles {{0,9,10},{10,18,19},{20,27,28},{0,12,13},{0,13,14},{0,14,15},{0,15,16},{0,1,16},{2,3,16},{1,2,16},{7,8,16},{6,7,16},{5,6,16},{4,5,16},{3,4,16},{9,10,17},{13,14,23},{12,13,23},{0,10,20},{10,19,20},{28,29,34},{0,11,22},{11,21,22},{0,12,22},{12,22,23},{15,16,23},{14,15,23},{10,17,24},{0,20,29},{0,11,29},{20,28,29},{18,19,25},{10,24,30},{19,20,26},{22,23,29},{11,21,29},{21,22,29},{20,27,37},{27,28,33},{28,34,38},{10,30,35},{19,25,31},{28,33,41},{20,26,32},{19,26,36},{20,37,43},{27,33,37},{20,32,43},{26,32,36},{19,31,36},{28,38,41},{37,41,43},{33,37,41},{10,35,39},{32,36,40},{32,40,43},{10,39,42},{10,42,44},{40,43,44},{36,40,44},{31,36,44},{25,31,44},{18,25,44},{10,18,44}};
  auto sign_triangles_pair = Utils::pcom_to_signs_and_triangles(pcom);
  REQUIRE(sign_triangles_pair.first == sign);
  REQUIRE(sign_triangles_pair.second == triangles);
  std::string pcom_created =  Utils::signs_and_triangles_to_pcom(sign_triangles_pair.first, sign_triangles_pair.second);
  auto sign_triangles_pair_created = Utils::pcom_to_signs_and_triangles(pcom_created);
  REQUIRE(sign_triangles_pair_created.first == sign);
  REQUIRE(sign_triangles_pair_created.second == triangles);
}

TEST_CASE("Isotopy::Graph from triangulation  (First odd degree Test)", "[isotopy_graph]") {
  int delta = 5;
  std::vector<bool> sign {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
  std::set<std::set<int>> triangles {{1, 0, 6}, {1, 6, 7}, {2, 1, 7}, {2, 7, 8}, {3, 2, 8}, {3, 8, 9}, {4, 3, 9}, {4, 9, 10}, {5, 4, 10}, {7, 6, 11}, {7, 11, 12}, {8, 7, 12}, {8, 12, 13}, {9, 8, 13}, {9, 13, 14}, {10, 9, 14}, {12, 11, 15}, {12, 15, 16}, {13, 12, 16}, {13, 16, 17}, {14, 13, 17}, {16, 15, 18}, {16, 18, 19}, {17, 16, 19}, {19, 18, 20}};
  Isotopy::Graph graph(delta, sign, triangles);
  graph.isotopy_type();
  std::string expected_viro = "<Jv1<1>>";
  REQUIRE(graph.viro_notation() == expected_viro);
  REQUIRE(graph.even_regions() == 2);
  REQUIRE(graph.odd_regions() == 1);
}

TEST_CASE("Isotopy::Graph from triangulation  (Test Case 0)", "[isotopy_graph]") {
    
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
    for (int i = 0; i < 10; ++i) {
      graph.isotopy_type();
    }
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


TEST_CASE("Isotopy::Graph 75", "[isotopy_graph]") {
    int delta = 8;
    std::vector<bool> sign {1,1,0,0,1,0,0,1,0,1,1,0,1,1,0,0,0,0,0,1,1,1,1,0,1,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0};
    std::set<std::set<int>> triangles {{0,9,10},{0,1,10},{9,10,17},{1,2,10},{10,17,18},{2,10,11},{10,18,19},{10,11,19},{17,18,25},{2,11,12},{18,19,25},{11,12,19},{17,24,25},{2,3,12},{24,25,30},{19,25,26},{12,19,20},{3,4,12},{19,26,27},{19,20,27},{25,26,32},{12,20,21},{25,30,31},{4,12,13},{25,31,32},{26,27,32},{20,21,27},{12,13,21},{30,31,36},{4,13,14},{31,32,36},{13,14,21},{27,32,33},{21,27,28},{30,35,36},{4,5,14},{27,33,34},{27,28,34},{32,36,37},{32,33,38},{21,28,29},{14,21,22},{35,36,39},{32,37,38},{21,22,29},{5,6,14},{33,34,38},{28,29,34},{36,37,41},{14,22,23},{37,38,41},{22,23,29},{36,39,40},{6,14,15},{36,40,41},{14,15,23},{39,40,43},{6,15,16},{40,41,43},{15,16,23},{39,42,43},{6,7,16},{42,43,44},{7,8,16}};
    Isotopy::Graph graph(delta, sign, triangles);
    graph.isotopy_type();
    REQUIRE(graph.even_regions() == 1);
    REQUIRE(graph.odd_regions() == 5);
    REQUIRE(graph.viro_notation() == "<1<5>>");
    //std::cout << "Shorthand: " << Utils::web_shorthand(delta, sign, triangles) << std::endl;
}

TEST_CASE("Graph constructor from edges (vector)", "[isotopy_graph]") {
  int delta = 2;
  std::vector<bool> sign = {true, true, false, false, false, true};
  std::vector<std::pair<int,int>> edges = {
      {0,1}, {1,2}, {0,3}, {3,4}, {2,4}, {1,3}, {1,4}, {4,5}, {3,5}
  };
  Isotopy::Graph graph(delta, sign, edges);
  graph.isotopy_type();
  REQUIRE(graph.viro_notation() == "<1>");
  REQUIRE(graph.even_regions() == 1);
  REQUIRE(graph.odd_regions() == 0);
}

TEST_CASE("Graph constructor from edges (set)", "[isotopy_graph]") {
  int delta = 2;
  std::vector<bool> sign = {true, true, false, false, false, true};
  std::set<std::pair<int,int>> edges = {
      {0,1}, {1,2}, {0,3}, {3,4}, {2,4}, {1,3}, {1,4}, {4,5}, {3,5}
  };
  Isotopy::Graph graph(delta, sign, edges);
  graph.isotopy_type();
  REQUIRE(graph.viro_notation() == "<1>");
  REQUIRE(graph.even_regions() == 1);
  REQUIRE(graph.odd_regions() == 0);
}

TEST_CASE("Component vs adjacency classification for delta=1", "[isotopy_graph]") {
  int delta = 1;
  std::vector<bool> sign {true, false, false};
  std::vector<Isotopy::Triangle> triangles = {{0,1,2}};

  Isotopy::Graph graph(delta, sign, triangles);
  REQUIRE(graph.component_edges.size() >= 3);
  REQUIRE(graph.adjacency_edges.size() >= 1);
}

TEST_CASE("Quadrant indices for delta=1", "[isotopy_graph]") {
  int delta = 1;
  std::vector<bool> sign {true, true, true};
  std::vector<Isotopy::Triangle> triangles = {{0,1,2}};

  Isotopy::Graph graph(delta, sign, triangles);

  REQUIRE(graph.quad_idxs[0][0] == 0);
  REQUIRE(graph.quad_idxs[0][1] == 0);
  REQUIRE(graph.quad_idxs[0][2] == 0);
  REQUIRE(graph.quad_idxs[0][3] == 0);

  REQUIRE(graph.quad_idxs[1][0] == 1);
  REQUIRE(graph.quad_idxs[1][1] == 1);
  REQUIRE(graph.quad_idxs[1][2] == 4);
  REQUIRE(graph.quad_idxs[1][3] == 4);

  REQUIRE(graph.quad_idxs[2][0] == 2);
  REQUIRE(graph.quad_idxs[2][1] == 3);
  REQUIRE(graph.quad_idxs[2][2] == 3);
  REQUIRE(graph.quad_idxs[2][3] == 2);

  REQUIRE(graph.antipodal_partner.size() == static_cast<std::size_t>(2 * (delta + 1)));
  auto vert_pair = graph.antipodal_partner[0];
  auto horiz_pair = graph.antipodal_partner[2];

  auto v0 = Isotopy::idx_to_point(delta, vert_pair.first);
  auto v1 = Isotopy::idx_to_point(delta, vert_pair.second);
  REQUIRE(v0.first == 0);
  REQUIRE(v0.second == 1);
  REQUIRE(v1.first == 0);
  REQUIRE(v1.second == -1);

  auto h0 = Isotopy::idx_to_point(delta, horiz_pair.first);
  auto h1 = Isotopy::idx_to_point(delta, horiz_pair.second);
  REQUIRE(h0.first == 1);
  REQUIRE(h0.second == 0);
  REQUIRE(h1.first == -1);
  REQUIRE(h1.second == 0);
}

TEST_CASE("Quadrant indices for delta=2", "[isotopy_graph]") {
  int delta = 2;
  std::vector<bool> sign {true, true, true, true, true, true};
  std::vector<Isotopy::Triangle> triangles = {{0,1,3}, {1,3,4}, {1,2,4}, {3,4,5}};

  Isotopy::Graph graph(delta, sign, triangles);

  REQUIRE(graph.quad_idxs[0][0] == 0);
  REQUIRE(graph.quad_idxs[0][1] == 0);
  REQUIRE(graph.quad_idxs[0][2] == 0);
  REQUIRE(graph.quad_idxs[0][3] == 0);

  REQUIRE(graph.quad_idxs[1][0] == 1);
  REQUIRE(graph.quad_idxs[1][1] == 1);
  REQUIRE(graph.quad_idxs[1][2] == 9);
  REQUIRE(graph.quad_idxs[1][3] == 9);

  REQUIRE(graph.quad_idxs[2][0] == 2);
  REQUIRE(graph.quad_idxs[2][1] == 2);
  REQUIRE(graph.quad_idxs[2][2] == 10);
  REQUIRE(graph.quad_idxs[2][3] == 10);

  REQUIRE(graph.quad_idxs[3][0] == 3);
  REQUIRE(graph.quad_idxs[3][1] == 6);
  REQUIRE(graph.quad_idxs[3][2] == 6);
  REQUIRE(graph.quad_idxs[3][3] == 3);

  REQUIRE(graph.quad_idxs[4][0] == 4);
  REQUIRE(graph.quad_idxs[4][1] == 7);
  REQUIRE(graph.quad_idxs[4][2] == 11);
  REQUIRE(graph.quad_idxs[4][3] == 12);

  REQUIRE(graph.quad_idxs[5][0] == 5);
  REQUIRE(graph.quad_idxs[5][1] == 8);
  REQUIRE(graph.quad_idxs[5][2] == 8);
  REQUIRE(graph.quad_idxs[5][3] == 5);

  REQUIRE(graph.antipodal_partner.size() == static_cast<std::size_t>(2 * (delta + 1)));
  auto pair_x0 = graph.antipodal_partner[0];
  auto pair_x1 = graph.antipodal_partner[2];
  auto pair_x2 = graph.antipodal_partner[4];

  auto x0a = Isotopy::idx_to_point(delta, pair_x0.first);
  auto x0b = Isotopy::idx_to_point(delta, pair_x0.second);
  REQUIRE(x0a.first == 0);
  REQUIRE(x0a.second == 2);
  REQUIRE(x0b.first == 0);
  REQUIRE(x0b.second == -2);

  auto x1a = Isotopy::idx_to_point(delta, pair_x1.first);
  auto x1b = Isotopy::idx_to_point(delta, pair_x1.second);
  REQUIRE(x1a.first == 1);
  REQUIRE(x1a.second == 1);
  REQUIRE(x1b.first == -1);
  REQUIRE(x1b.second == -1);

  auto x2a = Isotopy::idx_to_point(delta, pair_x2.first);
  auto x2b = Isotopy::idx_to_point(delta, pair_x2.second);
  REQUIRE(x2a.first == 2);
  REQUIRE(x2a.second == 0);
  REQUIRE(x2b.first == -2);
  REQUIRE(x2b.second == 0);
}


TEST_CASE("Canonical ordering - Curve1 (W&B case 1)", "[canonical_viro]") {
    int delta = 8;
    std::vector<bool> sign = {
        false, true, true, false, false, false, true, true, true, false, true, true, 
        true, true, true, true, true, true, false, true, false, true, false, true, 
        true, false, true, false, true, true, false, false, false, true, false, 
        true, true, true, false, false, false, false, false, true, false
    };
    std::set<std::set<int>> triangles = {
        {0,1,9}, {1,2,9}, {2,10,3}, {2,10,9}, {11,5,10}, {4,10,5}, {11,6,7}, {8,11,7}, 
        {3,10,4}, {11,5,6}, {12,8,13}, {8,13,14}, {22,14,15}, {11,19,12}, {8,22,16}, 
        {24,10,17}, {18,10,11}, {11,19,18}, {10,17,9}, {20,14,13}, {8,11,12}, {19,18,24}, 
        {20,14,21}, {14,15,8}, {12,20,13}, {20,21,24}, {8,22,15}, {32,16,28}, {24,10,18}, 
        {12,20,19}, {19,20,24}, {22,14,21}, {21,22,24}, {22,25,24}, {22,30,26}, {26,22,25}, 
        {25,30,24}, {29,39,33}, {16,33,23}, {27,30,31}, {26,30,25}, {31,30,35}, {35,16,22}, 
        {35,16,32}, {35,27,31}, {39,16,33}, {35,22,27}, {40,34,38}, {39,16,28}, {28,39,36}, 
        {22,30,27}, {36,28,32}, {36,32,35}, {29,39,37}, {37,29,34}, {33,29,23}, {36,39,35}, 
        {37,40,39}, {40,34,37}, {38,40,41}, {44,40,43}, {40,42,39}, {41,40,43}, {44,40,42}
    };
    
    Isotopy::Graph graph(delta, sign, triangles);
    graph.isotopy_type();
    
    // This should produce the canonical form <9v1<1>v1<10>>
    // (not the non-canonical <9v1<10>v1<1>>)
    REQUIRE(graph.even_regions() == 11);
    REQUIRE(graph.odd_regions() == 11);
    REQUIRE(graph.viro_notation() == "<9v1<1>v1<10>>");
}

TEST_CASE("Canonical ordering - Curve2 (W&B case 3)", "[canonical_viro]") {
    int delta = 8;
    std::vector<bool> sign = {
        false, true, false, true, false, true, false, true, true, false, true, true, 
        true, false, true, false, true, false, true, false, true, true, true, false, 
        false, false, false, false, true, true, false, true, false, false, false, 
        false, true, false, false, false, false, true, true, true, false
    };
    std::set<std::set<int>> triangles = {
        {1,17,9}, {2,11,3}, {2,10,11}, {1,17,10}, {5,11,4}, {2,10,1}, {5,12,11}, {6,18,12}, 
        {5,12,6}, {11,3,4}, {24,8,14}, {24,7,19}, {24,7,14}, {8,14,7}, {20,8,15}, {9,1,0}, 
        {10,11,18}, {18,11,12}, {10,17,18}, {6,18,13}, {19,18,24}, {6,13,7}, {19,7,13}, 
        {20,16,15}, {21,16,22}, {24,22,25}, {28,16,23}, {15,16,8}, {17,18,24}, {13,18,19}, 
        {24,8,20}, {24,16,20}, {24,22,21}, {16,35,32}, {24,16,21}, {26,22,25}, {25,31,30}, 
        {16,35,22}, {28,23,29}, {26,31,27}, {22,35,27}, {25,30,24}, {32,33,29}, {27,26,22}, 
        {32,29,28}, {33,43,41}, {31,35,30}, {42,43,44}, {41,33,38}, {32,16,28}, {36,37,39}, 
        {35,36,39}, {39,43,42}, {35,27,31}, {35,33,36}, {33,38,34}, {25,31,26}, {34,33,29}, 
        {36,37,33}, {37,39,40}, {43,37,40}, {35,33,32}, {33,43,37}, {39,43,40}
    };
    
    Isotopy::Graph graph(delta, sign, triangles);
    graph.isotopy_type();
    
    // This should produce the canonical form <5v1<1>v1<14>>
    // (not the non-canonical <5v1<14>v1<1>>)
    REQUIRE(graph.even_regions() == 7);
    REQUIRE(graph.odd_regions() == 15);
    REQUIRE(graph.viro_notation() == "<5v1<1>v1<14>>");
}

TEST_CASE("Isotopy::num_vertices - Basic calculation", "[isotopy]") {
    // Test the formula: nverts = (delta+1)(delta+2)/2
    REQUIRE(Isotopy::num_vertices(0) == 1);   // (0+1)(0+2)/2 = 1
    REQUIRE(Isotopy::num_vertices(1) == 3);   // (1+1)(1+2)/2 = 3
    REQUIRE(Isotopy::num_vertices(2) == 6);   // (2+1)(2+2)/2 = 6
    REQUIRE(Isotopy::num_vertices(3) == 10);  // (3+1)(3+2)/2 = 10
    REQUIRE(Isotopy::num_vertices(4) == 15);  // (4+1)(4+2)/2 = 15
    REQUIRE(Isotopy::num_vertices(5) == 21);  // (5+1)(5+2)/2 = 21
    REQUIRE(Isotopy::num_vertices(8) == 45);  // (8+1)(8+2)/2 = 45
}

TEST_CASE("Isotopy::num_edges - Basic calculation", "[isotopy]") {
    // Test the formula: nedges = 3*delta + 3*(delta^2 - delta)/2
    REQUIRE(Isotopy::num_edges(0) == 0);
    REQUIRE(Isotopy::num_edges(1) == 3);
    REQUIRE(Isotopy::num_edges(2) == 9);
    REQUIRE(Isotopy::num_edges(3) == 18);
    REQUIRE(Isotopy::num_edges(4) == 30);
    REQUIRE(Isotopy::num_edges(5) == 45);
}

TEST_CASE("Isotopy::triangles_to_edges - Basic triangle", "[isotopy]") {
    // Single triangle with vertices 0, 1, 2 for delta=2
    int delta = 2;
    std::set<std::set<int>> triangles = {{0, 1, 2}};

    auto edges = Isotopy::triangles_to_edges(delta, triangles);

    // Should produce 3 edges: (0,1), (0,2), (1,2)
    REQUIRE(edges.size() == 3);
    std::set<std::pair<int,int>> edge_set(edges.begin(), edges.end());
    REQUIRE(edge_set.count({0, 1}) == 1);
    REQUIRE(edge_set.count({0, 2}) == 1);
    REQUIRE(edge_set.count({1, 2}) == 1);
}

TEST_CASE("Isotopy::triangles_to_edges - Multiple triangles with shared edges", "[isotopy]") {
    // Two triangles sharing an edge: {0,1,2} and {1,2,3} for delta=3
    int delta = 3;
    std::set<std::set<int>> triangles = {{0, 1, 2}, {1, 2, 3}};

    auto edges = Isotopy::triangles_to_edges(delta, triangles);

    // Should produce 5 unique edges: (0,1), (0,2), (1,2), (1,3), (2,3)
    std::set<std::pair<int,int>> edge_set(edges.begin(), edges.end());
    REQUIRE(edge_set.size() == 5);
    REQUIRE(edge_set.count({0, 1}) == 1);
    REQUIRE(edge_set.count({0, 2}) == 1);
    REQUIRE(edge_set.count({1, 2}) == 1);
    REQUIRE(edge_set.count({1, 3}) == 1);
    REQUIRE(edge_set.count({2, 3}) == 1);
}

TEST_CASE("Isotopy::triangles_to_edges - Edge ordering consistency", "[isotopy]") {
    // Triangle with vertices in different orders for delta=8
    int delta = 8;
    std::set<std::set<int>> triangles = {{5, 2, 8}};

    auto edges = Isotopy::triangles_to_edges(delta, triangles);

    // Should produce edges with smaller index first
    std::set<std::pair<int,int>> edge_set(edges.begin(), edges.end());
    REQUIRE(edge_set.size() == 3);
    REQUIRE(edge_set.count({2, 5}) == 1);
    REQUIRE(edge_set.count({2, 8}) == 1);
    REQUIRE(edge_set.count({5, 8}) == 1);
}

TEST_CASE("Isotopy::triangles_to_edges - Invalid triangle", "[isotopy]") {
    // Triangle with wrong number of vertices
    int delta = 2;
    std::set<std::set<int>> triangles = {{0, 1, 2, 3}};

    REQUIRE_THROWS_AS(Isotopy::triangles_to_edges(delta, triangles), std::invalid_argument);
}
