#pragma once

#include <cstdlib>
#include <string>

// Data-dependent checks skip themselves when their input is absent, which makes
// a green run ambiguous: it may mean the check passed, or that it never ran.
// Setting the matching flag to 1 turns the skip into a failure, so a CI job that
// is supposed to have the data cannot pass by silently doing nothing.
//
//   ISOTOPY_REQUIRE_CORPUS   tests/isotopy_tests.yaml, shipped in this repo
//   ISOTOPY_REQUIRE_CATALOG  ISOTOPY_TYPE_INFO and ISOTOPY_REALIZABLE, from `searches`

inline std::string env_or_empty(const char* name) {
  const char* v = std::getenv(name);
  return v == nullptr ? std::string() : std::string(v);
}

inline bool env_flag(const char* name) { return env_or_empty(name) == "1"; }

#define ISOTOPY_SKIP_UNLESS_REQUIRED(flag, msg) \
  do {                                          \
    if (env_flag(flag)) FAIL(msg);              \
    WARN(msg);                                  \
  } while (0)
