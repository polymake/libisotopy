#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "isotopy_graph.h"
#include "catch.hpp"
#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <array>
#include <set>

#include <chrono>
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

TEST_CASE("Benchmark Isotopy::Graph batch test from file", "[benchmark][isotopy_graph]") {
    std::ifstream infile("tests/tree.txt");
    REQUIRE(infile);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(infile, line)) {
        lines.push_back(line);
    }
    REQUIRE(!lines.empty());

    int delta = 8;
    size_t sample_target = 10;
    std::vector<long long> durations;
    size_t case_count = 0;

    while (true) {
        for (const auto& l : lines) {
            auto data = parse_test_case(l);
            if (data.case_number == 1202044) {
                continue;
            }
            auto start = std::chrono::high_resolution_clock::now();
            Isotopy::Graph graph(delta, data.signs_vec, data.triangulation_vec);
            graph.isotopy_type();
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            durations.push_back(duration);
            ++case_count;
            if (durations.size() >= sample_target) break;
        }
    }

    // Analysis
    if (!durations.empty()) {
        // Mean
        double mean = std::accumulate(durations.begin(), durations.end(), 0.0) / durations.size();

        // Median
        std::vector<long long> sorted = durations;
        std::sort(sorted.begin(), sorted.end());
        double median;
        size_t n = sorted.size();
        if (n % 2 == 0)
            median = (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
        else
            median = sorted[n/2];

        // Stddev
        double sq_sum = std::inner_product(durations.begin(), durations.end(), durations.begin(), 0.0);
        double stddev = std::sqrt(sq_sum / n - mean * mean);

        std::cout << "Processed " << case_count << " cases." << std::endl;
        std::cout << "Mean: " << mean << " us, Median: " << median << " us, Stddev: " << stddev << " us." << std::endl;
    }
}
