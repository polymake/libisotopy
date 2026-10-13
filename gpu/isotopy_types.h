#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "isotopy_viro.h"

#ifdef ISOTOPY_HAVE_FKYAML
#include <fstream>
#include <iterator>

#include "node.hpp"
#endif

namespace IsotopyGPU {

struct TypeEntry {
  uint64_t code = 0ull;
  std::string viro;
  std::vector<std::string> aliases;
  int p = 0;
  int n = 0;
  int num_ovals = 0;
};

struct TypeTable {
  int delta = 0;
  std::vector<TypeEntry> entries;
  std::vector<uint64_t> keys;
  std::vector<std::string> rejected;
  uint64_t fingerprint = 0ull;

  uint16_t lookup(uint64_t code) const {
    const auto it = std::lower_bound(keys.begin(), keys.end(), code);
    if (it == keys.end() || *it != code) return kTypeUnknown;
    return static_cast<uint16_t>(it - keys.begin());
  }

  const TypeEntry* entry(uint16_t id) const {
    if (id >= entries.size()) return nullptr;
    return &entries[id];
  }
};

inline uint64_t mix64(uint64_t x) {
  x += 0x9E3779B97F4A7C15ull;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
  return x ^ (x >> 31);
}

inline TypeTable build_type_table(int delta, const std::vector<TypeEntry>& raw) {
  TypeTable table;
  table.delta = delta;
  const bool odd_degree = (delta % 2) != 0;

  for (const TypeEntry& in : raw) {
    TypeEntry out = in;
    try {
      const Viro::Parsed parsed = Viro::parse(in.viro);
      if (parsed.has_j != odd_degree) {
        table.rejected.push_back(in.viro);
        continue;
      }
      out.code = Viro::tree_code(parsed.tree);
      out.viro = Viro::to_viro(parsed.tree, parsed.has_j, false);
    } catch (const std::exception&) {
      table.rejected.push_back(in.viro);
      continue;
    }
    if (out.viro != in.viro) out.aliases.push_back(in.viro);
    table.entries.push_back(out);
  }

  std::sort(table.entries.begin(), table.entries.end(), [](const TypeEntry& a, const TypeEntry& b) {
    if (a.code != b.code) return a.code < b.code;
    return a.viro < b.viro;
  });

  std::vector<TypeEntry> merged;
  merged.reserve(table.entries.size());
  for (size_t i = 0; i < table.entries.size();) {
    size_t j = i;
    while (j < table.entries.size() && table.entries[j].code == table.entries[i].code) ++j;

    TypeEntry keep = table.entries[i];
    for (size_t k = i + 1; k < j; ++k) {
      const TypeEntry& other = table.entries[k];
      if (other.viro != keep.viro) {
        throw std::runtime_error("IsotopyGPU::build_type_table: degree " + std::to_string(delta) +
                                 " types \"" + keep.viro + "\" and \"" + other.viro +
                                 "\" share the tree code " + std::to_string(keep.code));
      }
      if (other.p != keep.p || other.n != keep.n || other.num_ovals != keep.num_ovals) {
        throw std::runtime_error("IsotopyGPU::build_type_table: degree " + std::to_string(delta) +
                                 " type \"" + keep.viro + "\" is listed twice with different " +
                                 "p/n/num_ovals");
      }
      keep.aliases.insert(keep.aliases.end(), other.aliases.begin(), other.aliases.end());
    }
    std::sort(keep.aliases.begin(), keep.aliases.end());
    keep.aliases.erase(std::unique(keep.aliases.begin(), keep.aliases.end()), keep.aliases.end());
    merged.push_back(keep);
    i = j;
  }
  table.entries = std::move(merged);

  if (table.entries.size() > static_cast<size_t>(kMaxTypeId) + 1) {
    throw std::runtime_error("IsotopyGPU::build_type_table: degree " + std::to_string(delta) +
                             " has " + std::to_string(table.entries.size()) +
                             " types, above the reserved-id ceiling");
  }

  table.keys.reserve(table.entries.size());
  for (const TypeEntry& e : table.entries) table.keys.push_back(e.code);

  uint64_t fingerprint = mix64(static_cast<uint64_t>(delta));
  for (uint64_t k : table.keys) fingerprint = mix64(fingerprint ^ k);
  table.fingerprint = fingerprint;
  return table;
}

#ifdef ISOTOPY_HAVE_FKYAML

inline std::map<int, std::vector<TypeEntry>> read_type_info(const std::string& yaml_path) {
  std::ifstream ifs(yaml_path);
  if (!ifs.is_open()) {
    throw std::runtime_error("IsotopyGPU::read_type_info: cannot open " + yaml_path);
  }
  const std::string text((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
  auto doc = fkyaml::node::deserialize(text);

  std::map<int, std::vector<TypeEntry>> raw;
  for (auto degree_it = doc.begin(); degree_it != doc.end(); ++degree_it) {
    const int delta = degree_it.key().get_value<int>();
    for (auto type_it = degree_it->begin(); type_it != degree_it->end(); ++type_it) {
      TypeEntry e;
      e.viro = type_it.key().get_value<std::string>();
      e.p = type_it->at("p").get_value<int>();
      e.n = type_it->at("n").get_value<int>();
      e.num_ovals = type_it->at("num_ovals").get_value<int>();
      raw[delta].push_back(e);
    }
  }
  return raw;
}

inline std::map<int, TypeTable> load_type_tables(const std::string& yaml_path) {
  const std::map<int, std::vector<TypeEntry>> raw = read_type_info(yaml_path);
  std::map<int, TypeTable> tables;
  for (const auto& kv : raw) tables[kv.first] = build_type_table(kv.first, kv.second);
  return tables;
}

inline TypeTable load_type_table(const std::string& yaml_path, int delta) {
  const std::map<int, std::vector<TypeEntry>> raw = read_type_info(yaml_path);
  const auto it = raw.find(delta);
  if (it == raw.end()) return build_type_table(delta, std::vector<TypeEntry>());
  return build_type_table(delta, it->second);
}

#endif

}  // namespace IsotopyGPU
