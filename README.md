# libisotopy

A C++ shared library that computes the **real scheme** (ambient isotopy type in
$\mathbb{RP}^2$) of a **T-curve** given by Viro's combinatorial patchworking: a
unimodular triangulation of $d\cdot\Delta_2$ together with a sign distribution
on its lattice points.

The algorithm runs in near-quadratic time in the degree $d$ via union-find on
the reflected signed triangulation, region-adjacency-graph traversal, and
rooted-tree canonicalization. It is described in:

> Geiselmann, Joswig, Kastner, Mundinger, Pokutta, Spiegel, Wack, Zimmer.
> *Fast Isotopy Computation for T-Curves.* ICMS 2026.
> [arXiv:2604.09221](https://arxiv.org/abs/2604.09221)

and was used to verify the realizability of all 121 degree-seven real schemes
([*121 Patchworked Curves of Degree Seven*](https://arxiv.org/abs/2602.06888);
data at [dmg-lab/CombinatorialPatchworking](https://github.com/dmg-lab/CombinatorialPatchworking)).

> **Status.** The core algorithm (`Isotopy::Graph`) is tested and stable. The
> surrounding build/packaging is geared toward our own use; if you embed the
> library elsewhere, please report any issues.

## Features

- `Isotopy::Graph` — computes the real scheme of a patchwork and emits its
  Viro (Rohlin–Viro) notation.
- Queries: `viro_notation()`, `even_regions()`, `odd_regions()`, `is_mcurve()`.

## Input format

`Isotopy::Graph(delta, sign, edges)` takes:
- `delta` — the degree $d$.
- `sign` — the sign distribution on the lattice points of $d\cdot\Delta_2$, in
  the standard lattice-point order (a `std::vector<bool>` of length
  $\tfrac12(d{+}1)(d{+}2)$).
- `edges` — the triangulation as an edge list over lattice-point indices
  (`std::vector<std::pair<int,int>>`).

## Documentation

API reference is provided as a man page:

```sh
man docs/man/man3/Isotopy_Graph.3
```

Worked examples are the test cases in `tests/isotopy_graph.cpp`.

## Building

A Makefile drives the build and tests:

```sh
make        # build the shared library
make test   # run tests
make clean  # remove build artifacts
```

## Usage

Include the header and link against the library:

```cpp
#include "isotopy_graph.h"
```
```sh
g++ myfile.cpp -L/path/to/lib -lisotopy
```

Example (degree 2):

```cpp
#include "isotopy_graph.h"
int main() {
    int delta = 2;
    std::vector<bool> sign = {true, true, true, true, true, true};
    std::vector<std::pair<int, int>> edges =
        {{0,1},{1,2},{0,3},{3,4},{2,4},{1,3},{1,4},{4,5},{3,5}};
    Isotopy::Graph g(delta, sign, edges);
    g.isotopy_type();
    std::string notation = g.viro_notation();
    int even = g.even_regions();
    int odd  = g.odd_regions();
    bool mcurve = g.is_mcurve();
}
```
