#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "isotopy_graph.h"

namespace Isotopy {

/// Lattice point (x, y) of the degree-delta simplex.
using Point = std::pair<int, int>;

/// Harnack split: 2 lattice points (simple) or 3 (double, apex in the middle).
using Split = std::vector<Point>;

/**
 * @brief Sign vector of the T-curve induced by a collection of Harnack splits.
 *
 * Harnack signs eta(x, y) = (x+1)(y+1) mod 2 twisted by each split: with
 * parity pair {a, b} the twist adds a0*b1+a1*b0 + (a1+b1)x + (a0+b0)y mod 2
 * on the split's zone, flipping exactly the two parity classes not in {a, b}.
 * Twists commute; pairwise compatibility is NOT checked.
 *
 * Zone convention: path canonicalized so its lexicographically smaller
 * endpoint comes first; zone = side with fewer simplex corners off the path;
 * tie goes to the side enclosed by the path plus the ccw boundary arc from
 * last endpoint back to first.  Result is independent of input vertex order.
 *
 * Split validity: endpoints on the boundary, parities alternate, edges
 * primitive, doubles non-collinear.
 *
 * @param delta The degree of the patchwork.
 * @param splits The split collection.
 * @return Signs of all lattice points in flat index order (matching Graph),
 *         length (delta+1)(delta+2)/2.
 * @throws std::invalid_argument on delta < 1 or an invalid split.
 */
std::vector<bool> signs_from_splits(int delta,
                                    const std::vector<Split>& splits);

/**
 * @brief Unimodular triangulation refining the splits.
 *
 * Specialization of refine(): the forced edges are the split edges.
 * Vertex indices match the sign vector of signs_from_splits.
 *
 * @param delta The degree of the patchwork.
 * @param splits The split collection.
 * @return The triangulation, or std::nullopt if the split edges conflict,
 *         e.g. for crossing splits.
 * @throws std::invalid_argument on delta < 1 or an invalid split.
 */
std::optional<std::vector<Triangle>> triangulation_from_splits(
    int delta, const std::vector<Split>& splits);

}
