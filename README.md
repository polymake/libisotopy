# libisotopy

A C++ shared library for computing isotopy types of patchworks.

## Documentation
For detailed API documentation, use the man page:  
```sh
man docs/man/man3/Isotopy_Graph.3
```

## Features

- Provides the `Isotopy::Graph` class for representing isotopy types of patchworks.
- Methods to compute isotopy types and generate Viro notation. 

## Examples

Examples can be found as test cases in the file `tests/isotopy_graph.cpp`.

## Building

This project uses a Makefile for building the shared library and running tests.

To build the shared library:

```sh
make
```

To run the tests:

```sh
make test
```

To clean build artifacts:

```sh
make clean
```

## Usage

1. Include the header in your project:
   ```cpp
   #include "isotopy_graph.h"
   ```

2. Link against the shared library:
   ```sh
   g++ myfile.cpp -L/path/to/lib -lisotopy
   ```

3. Example code:
```cpp
   #include "isotopy_graph.h"
   int main() {
       int delta = 2;
       std::vector<bool> sign = {true, true, true, true, true, true};
       std::set<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {0, 3}, {3, 4}, {2, 4}, {1,3}, {1,4}, {4,5}, {3,5}};
       Isotopy::Graph g(delta, sign, edges);
       g.isotopy_type();
       std::string notation = g.viro_notation();
       int even = g.even_regions();
       int odd = g.odd_regions(); 
       bool mcurve = g.is_mcurve();
   }
```

