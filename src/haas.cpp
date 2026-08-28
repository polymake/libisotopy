#include "haas.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>

namespace Isotopy {

namespace {

int cross(const Point& o, const Point& a, const Point& b) {
  return (a.first - o.first) * (b.second - o.second) -
         (a.second - o.second) * (b.first - o.first);
}

bool on_segment(const Point& a, const Point& b, const Point& q) {
  return std::min(a.first, b.first) <= q.first &&
         q.first <= std::max(a.first, b.first) &&
         std::min(a.second, b.second) <= q.second &&
         q.second <= std::max(a.second, b.second);
}

bool on_boundary(const Point& p, int delta) {
  return p.first == 0 || p.second == 0 || p.first + p.second == delta;
}

Point parity(const Point& p) { return {p.first % 2, p.second % 2}; }

// Row-major-in-y flat index, the parent-repo convention for Graph inputs.
// NOT point_to_idx, whose ordering is the transpose of this one; signs and
// triangulation must both use this convention to stay consistent.
int flat_idx(int delta, const Point& p) {
  return p.second * (delta + 1) - p.second * (p.second - 1) / 2 + p.first;
}

bool primitive(const Point& a, const Point& b) {
  return std::gcd(std::abs(a.first - b.first),
                  std::abs(a.second - b.second)) == 1;
}

std::string pt_str(const Point& p) {
  return "(" + std::to_string(p.first) + "," + std::to_string(p.second) + ")";
}

void validate_split(const Split& s, int delta) {
  if (s.size() != 2 && s.size() != 3)
    throw std::invalid_argument("split must have 2 or 3 vertices, got " +
                                std::to_string(s.size()));
  for (const Point& p : s)
    if (p.first < 0 || p.second < 0 || p.first + p.second > delta)
      throw std::invalid_argument("split vertex " + pt_str(p) +
                                  " outside the simplex");
  if (!on_boundary(s.front(), delta) || !on_boundary(s.back(), delta))
    throw std::invalid_argument("split endpoints must lie on the boundary");
  if (s.size() == 3) {
    if (parity(s[0]) != parity(s[2]) || parity(s[0]) == parity(s[1]))
      throw std::invalid_argument(
          "double split parities must satisfy p0 == p2 != p1");
    // A collinear double is a straight cut through its interior midpoint and
    // is a legal Harnack split; with the apex on the boundary the whole path
    // lies on the boundary and the zone degenerates.
    if (cross(s[0], s[2], s[1]) == 0 && on_boundary(s[1], delta))
      throw std::invalid_argument("collinear double split with boundary apex");
  } else if (parity(s[0]) == parity(s[1])) {
    throw std::invalid_argument("simple split endpoints have equal parity");
  }
  for (size_t k = 0; k + 1 < s.size(); ++k)
    if (!primitive(s[k], s[k + 1]))
      throw std::invalid_argument("non-primitive split edge " + pt_str(s[k]) +
                                  "-" + pt_str(s[k + 1]));
}

// All boundary lattice points of the simplex in ccw order.
std::vector<Point> boundary_walk(int delta) {
  std::vector<Point> walk;
  walk.reserve(3 * delta);
  for (int k = 0; k < delta; ++k) walk.push_back({k, 0});
  for (int k = 0; k < delta; ++k) walk.push_back({delta - k, k});
  for (int k = 0; k < delta; ++k) walk.push_back({0, delta - k});
  return walk;
}

// q inside or on the closed polygon; exact integer ray casting.
bool point_in_poly(const Point& q, const std::vector<Point>& poly) {
  const size_t n = poly.size();
  for (size_t k = 0; k < n; ++k) {
    const Point& a = poly[k];
    const Point& b = poly[(k + 1) % n];
    if (cross(a, b, q) == 0 && on_segment(a, b, q)) return true;
  }
  bool inside = false;
  for (size_t k = 0; k < n; ++k) {
    const Point& a = poly[k];
    const Point& b = poly[(k + 1) % n];
    if ((a.second > q.second) != (b.second > q.second)) {
      // x of the edge at height q.y exceeds q.x iff num/dy > 0
      const int dy = b.second - a.second;
      const int num = (a.first - q.first) * dy +
                      (q.second - a.second) * (b.first - a.first);
      if (num * dy > 0) inside = !inside;
    }
  }
  return inside;
}

std::vector<Point> lattice_in_poly(const std::vector<Point>& poly, int delta) {
  std::vector<Point> out;
  for (int y = 0; y <= delta; ++y)
    for (int x = 0; x <= delta - y; ++x)
      if (point_in_poly({x, y}, poly)) out.push_back({x, y});
  return out;
}

// Zone of a canonically oriented split, per the convention in haas.h.
std::vector<Point> zone_of(const Split& s, int delta) {
  const std::vector<Point> walk = boundary_walk(delta);
  const auto idx_of = [&](const Point& p) {
    return static_cast<size_t>(std::find(walk.begin(), walk.end(), p) -
                               walk.begin());
  };
  const size_t i0 = idx_of(s.front());
  const size_t i1 = idx_of(s.back());
  const auto arc = [&](size_t a, size_t b) {
    std::vector<Point> out;
    for (size_t k = a; k != b; k = (k + 1) % walk.size()) out.push_back(walk[k]);
    out.push_back(walk[b]);
    return out;
  };
  const std::vector<Point> arc1 = arc(i1, i0);
  const std::vector<Point> arc2 = arc(i0, i1);
  // the path closed with either boundary arc gives the two candidate zones
  std::vector<Point> poly1(s.begin(), s.end() - 1);
  poly1.insert(poly1.end(), arc1.begin(), arc1.end() - 1);
  std::vector<Point> poly2(s.rbegin(), s.rend() - 1);
  poly2.insert(poly2.end(), arc2.begin(), arc2.end() - 1);
  const std::vector<Point> z1 = lattice_in_poly(poly1, delta);
  const std::vector<Point> z2 = lattice_in_poly(poly2, delta);
  const std::vector<Point> corners = {{0, 0}, {delta, 0}, {0, delta}};
  const auto corners_off_path = [&](const std::vector<Point>& z) {
    int c = 0;
    for (const Point& p : corners)
      if (std::find(z.begin(), z.end(), p) != z.end() &&
          std::find(s.begin(), s.end(), p) == s.end())
        ++c;
    return c;
  };
  return corners_off_path(z1) <= corners_off_path(z2) ? z1 : z2;
}

}  // namespace

std::vector<bool> signs_from_splits(int delta,
                                    const std::vector<Split>& splits) {
  if (delta < 1) throw std::invalid_argument("delta must be >= 1");
  std::vector<bool> sigma((delta + 1) * (delta + 2) / 2);
  for (int y = 0; y <= delta; ++y)
    for (int x = 0; x <= delta - y; ++x)
      sigma[flat_idx(delta, {x, y})] = ((x + 1) * (y + 1)) % 2;
  for (Split s : splits) {
    validate_split(s, delta);
    if (s.back() < s.front()) std::reverse(s.begin(), s.end());
    const Point a = s.size() == 3 ? parity(s[1]) : parity(s[0]);
    const Point b = s.size() == 3 ? parity(s[0]) : parity(s[1]);
    const int eps = (a.first * b.second + a.second * b.first) % 2;
    const int i = (a.second + b.second) % 2;
    const int j = (a.first + b.first) % 2;
    for (const Point& q : zone_of(s, delta)) {
      const int flip = (eps + i * q.first + j * q.second) % 2;
      if (flip) sigma[flat_idx(delta, q)] = !sigma[flat_idx(delta, q)];
    }
  }
  return sigma;
}

bool compatible(int delta, const std::vector<Split>& splits) {
  if (delta < 1) throw std::invalid_argument("delta must be >= 1");
  std::vector<std::pair<Point, Point>> edges;
  for (const Split& s : splits) {
    validate_split(s, delta);
    for (size_t k = 0; k + 1 < s.size(); ++k)
      edges.push_back(std::minmax(s[k], s[k + 1]));
  }
  for (size_t a = 0; a < edges.size(); ++a)
    for (size_t b = a + 1; b < edges.size(); ++b) {
      const auto& [p, q] = edges[a];
      const auto& [r, t] = edges[b];
      if (edges[a] == edges[b]) return false;
      const long s1 = (long)cross(p, q, r) * cross(p, q, t);
      const long s2 = (long)cross(r, t, p) * cross(r, t, q);
      if (s1 < 0 && s2 < 0) return false;
    }
  return true;
}

std::optional<std::vector<Triangle>> triangulation_from_splits(
    int delta, const std::vector<Split>& splits) {
  if (!compatible(delta, splits)) return std::nullopt;
  std::vector<Edge> edges;
  for (const Split& s : splits)
    for (size_t k = 0; k + 1 < s.size(); ++k)
      edges.push_back({flat_idx(delta, s[k]), flat_idx(delta, s[k + 1])});
  return refine(delta, edges);
}

}
