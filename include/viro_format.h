#pragma once

#include <string>
#include <vector>

namespace ViroFormat {

using std::vector;
using std::string;

struct NodeSummary {
  int leaf_count = 0;
  vector<const string*> children;  // pointers remain valid for the call duration
};

string format(const NodeSummary& summary, bool unicode = false);

}  // namespace ViroFormat

