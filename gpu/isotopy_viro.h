#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "isotopy_gpu.cuh"

namespace IsotopyGPU {
namespace Viro {

inline constexpr const char* kOpenAscii = "<";
inline constexpr const char* kCloseAscii = ">";
inline constexpr const char* kSepAscii = "v";
inline constexpr const char* kOpenUnicode = "⟨";
inline constexpr const char* kCloseUnicode = "⟩";
inline constexpr const char* kSepUnicode = "⊔";
inline constexpr const char* kEmptySet = "∅";
inline constexpr int kMaxLeafCount = 1 << 16;

struct Tree {
  int leaf_count = 0;
  std::vector<Tree> children;
};

struct Parsed {
  Tree tree;
  bool has_j = false;
};

namespace detail {

class Reader {
 public:
  explicit Reader(const std::string& text) : s_(text) {}

  Parsed run() {
    Parsed out;
    if (s_ == kEmptySet) return out;
    expect_open();
    if (peek_byte('J')) {
      ++i_;
      out.has_j = true;
      match_sep();
    }
    out.tree = inner();
    expect_close();
    if (i_ != s_.size()) fail("trailing characters");
    return out;
  }

 private:
  [[noreturn]] void fail(const std::string& what) const {
    throw std::invalid_argument("IsotopyGPU::Viro::parse: " + what + " at position " +
                                std::to_string(i_) + " of \"" + s_ + "\"");
  }

  bool match(const char* token) {
    const size_t n = std::string(token).size();
    if (s_.compare(i_, n, token) != 0) return false;
    i_ += n;
    return true;
  }

  bool looks_at(const char* token) const { return s_.compare(i_, std::string(token).size(), token) == 0; }

  bool peek_byte(char c) const { return i_ < s_.size() && s_[i_] == c; }

  bool match_open() { return match(kOpenAscii) || match(kOpenUnicode); }
  bool match_close() { return match(kCloseAscii) || match(kCloseUnicode); }
  bool match_sep() { return match(kSepAscii) || match(kSepUnicode); }
  bool looks_at_open() const { return looks_at(kOpenAscii) || looks_at(kOpenUnicode); }
  bool looks_at_close() const { return looks_at(kCloseAscii) || looks_at(kCloseUnicode); }

  void expect_open() {
    if (!match_open()) fail("expected an opening delimiter");
  }
  void expect_close() {
    if (!match_close()) fail("expected a closing delimiter");
  }

  int read_int() {
    const size_t start = i_;
    long value = 0;
    while (i_ < s_.size() && s_[i_] >= '0' && s_[i_] <= '9') {
      value = value * 10 + (s_[i_] - '0');
      if (value > kMaxLeafCount) fail("count exceeds " + std::to_string(kMaxLeafCount));
      ++i_;
    }
    if (i_ == start) fail("expected a count");
    return static_cast<int>(value);
  }

  void read_group(Tree& parent, int repeat) {
    if (repeat <= 0) fail("group repeat count must be positive");
    expect_open();
    Tree sub = inner();
    expect_close();
    for (int k = 0; k < repeat; ++k) parent.children.push_back(sub);
  }

  Tree inner() {
    Tree t;
    if (i_ >= s_.size() || looks_at_close()) return t;

    const int first = read_int();
    if (looks_at_open()) {
      read_group(t, first);
    } else {
      t.leaf_count = first;
      if (!match_sep()) return t;
      const int repeat = read_int();
      read_group(t, repeat);
    }
    while (match_sep()) {
      const int repeat = read_int();
      read_group(t, repeat);
    }
    return t;
  }

  const std::string& s_;
  size_t i_ = 0;
};

struct Coded {
  uint64_t code;
  int len;
};

inline bool coded_less(const Coded& a, const Coded& b) {
  if (a.len != b.len) return a.len < b.len;
  return a.code < b.code;
}

inline Coded code_of(const Tree& t) {
  std::vector<Coded> kids;
  kids.reserve(static_cast<size_t>(t.leaf_count) + t.children.size());
  for (int k = 0; k < t.leaf_count; ++k) kids.push_back(Coded{0ull, 0});
  for (const Tree& c : t.children) kids.push_back(code_of(c));
  std::sort(kids.begin(), kids.end(), coded_less);

  uint64_t code = 0ull;
  int off = 0;
  for (const Coded& k : kids) {
    if (off + k.len + 2 > kCodeBits) throw std::overflow_error("IsotopyGPU::Viro::tree_code");
    code |= (1ull | (k.code << 1)) << off;
    off += k.len + 2;
  }
  return Coded{code, off};
}

inline void count_by_depth(const Tree& t, int depth, int* even, int* odd) {
  ++(*(depth % 2 == 0 ? even : odd));
  for (int k = 0; k < t.leaf_count; ++k) ++(*((depth + 1) % 2 == 0 ? even : odd));
  for (const Tree& c : t.children) count_by_depth(c, depth + 1, even, odd);
}

inline bool is_leaf(const Tree& t) { return t.leaf_count == 0 && t.children.empty(); }

inline std::string render(const Tree& t, const char* open, const char* close, const char* sep) {
  int leaf_count = t.leaf_count;
  std::vector<const Tree*> non_leaf;
  for (const Tree& c : t.children) {
    if (is_leaf(c)) {
      ++leaf_count;
    } else {
      non_leaf.push_back(&c);
    }
  }

  if (non_leaf.empty()) {
    if (leaf_count == 0) return std::string(open) + close;
    return std::string(open) + std::to_string(leaf_count) + close;
  }
  if (leaf_count == 0 && non_leaf.size() == 1) {
    return std::string(open) + "1" + render(*non_leaf.front(), open, close, sep) + close;
  }

  std::map<std::string, int> counts;
  for (const Tree* c : non_leaf) ++counts[render(*c, open, close, sep)];

  std::vector<std::pair<int, std::string>> grouped;
  grouped.reserve(counts.size());
  for (const auto& kv : counts) grouped.push_back(std::make_pair(kv.second, kv.first));
  std::sort(grouped.begin(), grouped.end(),
            [](const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
              if (a.second.length() != b.second.length()) return a.second.length() < b.second.length();
              return a.second < b.second;
            });

  std::string result = open;
  if (leaf_count > 0) result += std::to_string(leaf_count) + sep;
  for (size_t g = 0; g < grouped.size(); ++g) {
    if (g != 0) result += sep;
    result += std::to_string(grouped[g].first) + grouped[g].second;
  }
  result += close;
  return result;
}

}  // namespace detail

inline Parsed parse(const std::string& s) { return detail::Reader(s).run(); }

inline uint64_t tree_code(const Tree& t) {
  const detail::Coded root = detail::code_of(t);
  if (root.len >= kCodeBits) throw std::overflow_error("IsotopyGPU::Viro::tree_code");
  return (1ull << root.len) | root.code;
}

inline int node_count(const Tree& t) {
  int n = 1 + t.leaf_count;
  for (const Tree& c : t.children) n += node_count(c);
  return n;
}

inline int num_ovals(const Tree& t, bool has_j) { return node_count(t) - 1 + (has_j ? 1 : 0); }

inline void region_pn(const Tree& t, bool odd_degree, int* p, int* n) {
  int even = 0;
  int odd = 0;
  detail::count_by_depth(t, 0, &even, &odd);
  if (odd_degree) {
    *p = even;
    *n = odd;
  } else {
    *p = odd;
    *n = even - 1;
  }
}

inline std::string to_viro(const Tree& t, bool has_j, bool unicode) {
  const char* open = unicode ? kOpenUnicode : kOpenAscii;
  const char* close = unicode ? kCloseUnicode : kCloseAscii;
  const char* sep = unicode ? kSepUnicode : kSepAscii;

  std::string notation = detail::render(t, open, close, sep);
  if (!has_j) return notation;

  const size_t open_len = std::string(open).size();
  const bool has_additional_components = notation.size() > open_len + std::string(close).size();
  notation.insert(open_len, has_additional_components ? std::string("J") + sep : std::string("J"));
  return notation;
}

}  // namespace Viro
}  // namespace IsotopyGPU
