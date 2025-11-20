#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include "node.hpp"
#include "isotopy_graph.h"

#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <functional>

TEST_CASE("Benchmark Sebastian's example", "[isotopy_graph]") {
  int delta = 6;
  std::vector<bool> sign {0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,0,0,1,0,1,0,1,1,1,1,0,1,1};
  std::set<std::set<int>> triangles {{5,6,12},{5,11,12},{4,5,11},{11,12,17},{4,10,11},{11,16,17},{3,4,10},{10,11,16},{16,17,21},{2,3,10},{2,9,10},{9,10,16},{1,2,9},{9,16,21},{1,9,21},{1,15,21},{1,8,15},{1,7,8},{7,8,15},{0,1,7},{7,15,21},{7,14,21},{14,20,21},{14,19,20},{7,13,14},{13,14,19},{20,21,24},{20,23,24},{19,20,23},{23,24,26},{13,18,19},{19,22,23},{18,19,22},{23,25,26},{22,23,25},{25,26,27}};
  BENCHMARK_ADVANCED("constructor")(Catch::Benchmark::Chronometer meter) {
    meter.measure([&] { return Isotopy::Graph(delta, sign, triangles); });
  };
  BENCHMARK_ADVANCED("isotopy_type")(Catch::Benchmark::Chronometer meter) {
    Isotopy::Graph graph(delta, sign, triangles);
    meter.measure([&] { return graph.isotopy_type(); });
  };
  BENCHMARK_ADVANCED("viro_notation")(Catch::Benchmark::Chronometer meter) {
    Isotopy::Graph graph(delta, sign, triangles);
    graph.isotopy_type();
    meter.measure([&] { return graph.viro_notation(); });
  };
}

TEST_CASE("Benchmark specific cases from YAML (triangles)", "[isotopy_graph][yaml][specific][triangles]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  std::vector<int> case_numbers = {21, 1358, 1351};
  for (int case_num : case_numbers) {
    if (case_num - 1 < 0 || static_cast<size_t>(case_num - 1) >= doc.size()) continue;
    auto& node = doc[case_num - 1];
    int delta = node.at("degree").get_value<int>();
    std::vector<bool> signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
    std::set<std::set<int>> triangulation_vec = node.at("triangulation").get_value<std::set<std::set<int>>>();

    BENCHMARK_ADVANCED("triangles: constructor " + std::to_string(case_num))(Catch::Benchmark::Chronometer meter) {
      meter.measure([&] { return Isotopy::Graph(delta, signs_vec, triangulation_vec); });
    };


    BENCHMARK_ADVANCED("triangles: isotopy_type " + std::to_string(case_num))(Catch::Benchmark::Chronometer meter) {
      Isotopy::Graph graph(delta, signs_vec, triangulation_vec);
      meter.measure([&] { return graph.isotopy_type(); });
    };

    BENCHMARK_ADVANCED("triangles: viro_notation " + std::to_string(case_num))(Catch::Benchmark::Chronometer meter) {
      Isotopy::Graph graph(delta, signs_vec, triangulation_vec);
      graph.isotopy_type();
      meter.measure([&] { return graph.viro_notation(); });
    };
    BENCHMARK_ADVANCED("triangles: all three " + std::to_string(case_num))(Catch::Benchmark::Chronometer meter) {
      meter.measure([&] {
        Isotopy::Graph graph(delta, signs_vec, triangulation_vec);
        return graph.viro_notation();
      });
    };
  }
}

TEST_CASE("Benchmark specific cases from YAML (edges)", "[isotopy_graph][yaml][specific][edges]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  std::vector<int> case_numbers = {21, 1358, 1351};
  for (int case_num : case_numbers) {
    if (case_num - 1 < 0 || static_cast<size_t>(case_num - 1) >= doc.size()) continue;
    auto& node = doc[case_num - 1];
    int delta = node.at("degree").get_value<int>();
    std::vector<bool> signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
    std::vector<std::pair<int, int>> edges_vec = node.at("edges").get_value<std::vector<std::pair<int, int>>>();

    BENCHMARK_ADVANCED("edges: constructor " + std::to_string(case_num))(Catch::Benchmark::Chronometer meter) {
      meter.measure([&] { return Isotopy::Graph(delta, signs_vec, edges_vec); });
    };

    BENCHMARK_ADVANCED("edges: isotopy_type " + std::to_string(case_num))(Catch::Benchmark::Chronometer meter) {
      Isotopy::Graph graph(delta, signs_vec, edges_vec);
      meter.measure([&] { return graph.isotopy_type(); });
    };

    BENCHMARK_ADVANCED("edges: viro_notation " + std::to_string(case_num))(Catch::Benchmark::Chronometer meter) {
      Isotopy::Graph graph(delta, signs_vec, edges_vec);
      graph.isotopy_type();
      meter.measure([&] { return graph.viro_notation(); });
    };

    BENCHMARK_ADVANCED("edges: all three " + std::to_string(case_num))(Catch::Benchmark::Chronometer meter) {
      meter.measure([&] {
        Isotopy::Graph graph(delta, signs_vec, edges_vec);
        return graph.viro_notation();
      });
    };
  }
}

TEST_CASE("Benchmark single case from YAML for 10 seconds (triangles, averaged over 3 runs)", "[isotopy_graph][yaml][timing][triangles]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  int case_num = 1358;
  auto& node = doc[case_num - 1];
  int delta = node.at("degree").get_value<int>();
  std::vector<bool> signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
  std::set<std::set<int>> triangulation_vec = node.at("triangulation").get_value<std::set<std::set<int>>>();

  using clock = std::chrono::steady_clock;
  int total_count = 0;
  int runs = 3;
  for (int i = 0; i < runs; ++i) {
    auto start = clock::now();
    int count = 0;
    while (std::chrono::duration_cast<std::chrono::seconds>(clock::now() - start).count() < 10) {
      Isotopy::Graph graph(delta, signs_vec, triangulation_vec);
      volatile std::string viro = graph.viro_notation();
      ++count;
    }
    std::cout << "Run " << (i + 1) << ": Computed isotopy_type " << count << " times in 10 seconds.\n";
    total_count += count;
  }
  double average = static_cast<double>(total_count) / runs;
  std::cout << "Average: Computed isotopy_type " << average << " times in 10 seconds (over " << runs << " runs).\n";
}

TEST_CASE("Benchmark single case from YAML for 10 seconds (edges, averaged over 3 runs)", "[isotopy_graph][yaml][timing][edges]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  int case_num = 1358;
  auto& node = doc[case_num - 1];
  int delta = node.at("degree").get_value<int>();
  std::vector<bool> signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
  std::vector<std::pair<int, int>> edges_vec = node.at("edges").get_value<std::vector<std::pair<int, int>>>();

  using clock = std::chrono::steady_clock;
  int total_count = 0;
  int runs = 3;
  for (int i = 0; i < runs; ++i) {
    auto start = clock::now();
    int count = 0;
    while (std::chrono::duration_cast<std::chrono::seconds>(clock::now() - start).count() < 10) {
      Isotopy::Graph graph(delta, signs_vec, edges_vec);
      volatile std::string viro = graph.viro_notation();
      ++count;
    }
    std::cout << "Run " << (i + 1) << ": Computed isotopy_type " << count << " times in 10 seconds.\n";
    total_count += count;
  }
  double average = static_cast<double>(total_count) / runs;
  std::cout << "Average: Computed isotopy_type " << average << " times in 10 seconds (over " << runs << " runs).\n";
}

TEST_CASE("Benchmark sparse case from YAML for 10 seconds (triangles, averaged over 3 runs)", "[isotopy_graph][yaml][timing][triangles][sparse]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  int case_num = 21;
  if (case_num - 1 < 0 || static_cast<size_t>(case_num - 1) >= doc.size()) {
    FAIL("Sparse benchmark case not found in YAML");
  }
  auto& node = doc[case_num - 1];
  int delta = node.at("degree").get_value<int>();
  std::vector<bool> signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
  std::set<std::set<int>> triangulation_vec = node.at("triangulation").get_value<std::set<std::set<int>>>();

  using clock = std::chrono::steady_clock;
  int total_count = 0;
  int runs = 3;
  for (int i = 0; i < runs; ++i) {
    auto start = clock::now();
    int count = 0;
    while (std::chrono::duration_cast<std::chrono::seconds>(clock::now() - start).count() < 10) {
      Isotopy::Graph graph(delta, signs_vec, triangulation_vec);
      volatile std::string viro = graph.viro_notation();
      ++count;
    }
    std::cout << "Run " << (i + 1) << ": Computed isotopy_type " << count << " times in 10 seconds.\n";
    total_count += count;
  }
  double average = static_cast<double>(total_count) / runs;
  std::cout << "Average: Computed isotopy_type " << average << " times in 10 seconds (over " << runs << " runs).\n";
}

TEST_CASE("Benchmark sparse case from YAML for 10 seconds (edges, averaged over 3 runs)", "[isotopy_graph][yaml][timing][edges][sparse]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  int case_num = 21;
  if (case_num - 1 < 0 || static_cast<size_t>(case_num - 1) >= doc.size()) {
    FAIL("Sparse benchmark case not found in YAML");
  }
  auto& node = doc[case_num - 1];
  int delta = node.at("degree").get_value<int>();
  std::vector<bool> signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
  std::vector<std::pair<int, int>> edges_vec = node.at("edges").get_value<std::vector<std::pair<int, int>>>();

  using clock = std::chrono::steady_clock;
  int total_count = 0;
  int runs = 3;
  for (int i = 0; i < runs; ++i) {
    auto start = clock::now();
    int count = 0;
    while (std::chrono::duration_cast<std::chrono::seconds>(clock::now() - start).count() < 10) {
      Isotopy::Graph graph(delta, signs_vec, edges_vec);
      volatile std::string viro = graph.viro_notation();
      ++count;
    }
    std::cout << "Run " << (i + 1) << ": Computed isotopy_type " << count << " times in 10 seconds.\n";
    total_count += count;
  }
  double average = static_cast<double>(total_count) / runs;
  std::cout << "Average: Computed isotopy_type " << average << " times in 10 seconds (over " << runs << " runs).\n";
}

TEST_CASE("Benchmark PCOM conversions", "[utils][pcom][benchmark]") {
  static const std::string sample_pcom = R"({
    "_ns": { "polymake": [ "https://polymake.org", "4.13" ] },
    "_type": "tropical::Hypersurface<Min>",
    "_id": "benchmark",
    "_libisotopy_version": "2",
    "_attrs": { "TYPE": { "attachment": true } },
    "TYPE": "<13v1<6>v1<1>>",
    "DUAL_SUBDIVISION": {
      "MAXIMAL_CELLS": [
        [0,9,10],[10,18,19],[20,27,28],[0,12,13],[0,13,14],[0,14,15],[0,15,16],[0,1,16],
        [2,3,16],[1,2,16],[7,8,16],[6,7,16],[5,6,16],[4,5,16],[3,4,16],[9,10,17],[13,14,23],
        [12,13,23],[0,10,20],[10,19,20],[28,29,34],[0,11,22],[11,21,22],[0,12,22],[12,22,23],
        [15,16,23],[14,15,23],[10,17,24],[0,20,29],[0,11,29],[20,28,29],[18,19,25],[10,24,30],
        [19,20,26],[22,23,29],[11,21,29],[21,22,29],[20,27,37],[27,28,33],[28,29,38],[10,30,35],
        [19,25,31],[28,33,41],[20,26,32],[19,26,36],[20,37,43],[27,33,37],[20,32,43],[26,32,36],
        [19,31,36],[28,38,41],[37,41,43],[33,37,41],[10,35,39],[32,36,40],[32,40,43],[10,39,42],
        [10,42,44],[40,43,44],[36,40,44],[31,36,44],[25,31,44],[18,25,44],[10,18,44]
      ],
      "WEIGHTS": []
    },
    "PATCHWORK": [{
      "_id": "benchmark#0",
      "SIGNS": [
        true,true,true,true,true,true,true,true,false,true,false,true,true,true,true,true,
        true,true,true,true,true,true,true,false,true,false,true,true,true,true,true,true,true,true,false,true,true,true,true,true,false,false,true,true
      ]
    }]
  })";

  const auto baseline = Utils::pcom_to_signs_and_triangles(sample_pcom);

  BENCHMARK("pcom_to_signs_and_triangles") {
    auto parsed = Utils::pcom_to_signs_and_triangles(sample_pcom);
    return parsed.first.size() + parsed.second.size();
  };

  BENCHMARK("signs_and_triangles_to_pcom") {
    return Utils::signs_and_triangles_to_pcom(baseline.first, baseline.second, "benchmark").size();
  };
}


/*
TEST_CASE("Isotopy::Graph batch timing analysis from YAML file", "[isotopy_graph][yaml][timing]") {
  std::ifstream ifs("tests/isotopy_tests.yaml");
  std::string yaml((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(yaml);

  // Store timings and case numbers
  std::vector<std::pair<double, int>> timings(doc.size(), {0.0, 0});

  int repeat_count = 50;
  for (int rep = 0; rep < repeat_count; ++rep) {
    int line_count = 0;
    if (rep % 5 == 0) {
      std::cout << "Starting repetition " << (rep + 1) << " of " << repeat_count << "\n";
    }
    for (auto& node : doc) {
      ++line_count;

      int delta = node.at("degree").get_value<int>();
      std::vector<bool> signs_vec = node.at("polarisation").get_value<std::vector<bool>>();
      std::set<std::set<int>> triangulation_vec = node.at("triangulation").get_value<std::set<std::set<int>>>();
      size_t nverts = Isotopy::num_vertices(delta);
      if (signs_vec.size() != nverts) {
        continue; // Skip invalid test case
      }

      auto start = std::chrono::high_resolution_clock::now();
      Isotopy::Graph graph(delta, signs_vec, triangulation_vec);
      graph.isotopy_type();
      graph.viro_notation();
      auto end = std::chrono::high_resolution_clock::now();
      double elapsed = std::chrono::duration<double, std::milli>(end - start).count();

      // Accumulate timings by index
      int idx = line_count - 1;
      timings[idx].first += elapsed;
      timings[idx].second = line_count;
    }
  }

  // Compute average
  for (auto& t : timings) {
    t.first /= repeat_count;
  }

  // Sort timings descending and print the worst 10
std::sort(timings.begin(), timings.end(), std::greater<std::pair<double, int>>());
  std::cout << "Worst 10 cases (elapsed ms, case number):\n";
  for (size_t i = 0; i < std::min<size_t>(10, timings.size()); ++i) {
    std::cout << timings[i].first << " ms, case #" << timings[i].second << "\n";
  }

}
*/



/*
0.283714 ms, case #1551
0.283089 ms, case #1358
0.282703 ms, case #1393
0.281751 ms, case #1443
0.281654 ms, case #1395
0.281499 ms, case #1400
0.281202 ms, case #1353
0.280646 ms, case #1308
0.2799 ms, case #1351
0.279849 ms, case #1397
*/
