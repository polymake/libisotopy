#include "viro_format.h"

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

using std::string;
using std::pair;
using std::map;
using std::vector;

struct Symbols {
  const char* open;
  const char* close;
  const char* separator;
};

Symbols symbols(bool unicode) {
  if (unicode) {
    return {"\u27E8", "\u27E9", "\u2294"};
  }
  return {"<", ">", "v"};
}

}  // namespace

namespace ViroFormat {

string format(const NodeSummary& summary, bool unicode) {
  const auto sym = symbols(unicode);

  if (summary.leaf_count == 0 && summary.children.empty()) {
    return string(sym.open) + sym.close;
  }

  if (summary.children.empty()) {
    return string(sym.open) + std::to_string(summary.leaf_count) + sym.close;
  }

  if (summary.leaf_count == 0 && summary.children.size() == 1) {
    return string(sym.open) + "1" + *summary.children.front() + sym.close;
  }

  map<string, int> grouped;
  for (const auto* child : summary.children) {
    grouped[*child]++;
  }

  vector<pair<int, string>> ordered;
  ordered.reserve(grouped.size());
  for (const auto& entry : grouped) {
    ordered.emplace_back(entry.second, entry.first);
  }

  std::sort(ordered.begin(), ordered.end(),
            [](const auto& a, const auto& b) {
              if (a.second.size() != b.second.size()) {
                return a.second.size() < b.second.size();
              }
              return a.second < b.second;
            });

  string result(sym.open);
  if (summary.leaf_count > 0) {
    result += std::to_string(summary.leaf_count);
    result += sym.separator;
  }
  for (std::size_t idx = 0; idx < ordered.size(); ++idx) {
    if (idx > 0) {
      result += sym.separator;
    }
    result += std::to_string(ordered[idx].first);
    result += ordered[idx].second;
  }
  result += sym.close;
  return result;
}

}  // namespace ViroFormat
