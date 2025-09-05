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
       Isotopy::Graph g(3, {true, false, true}, {{0,1},{1,2}});
       g.isotopy_type();
       std::string notation = g.viro_notation();
       int even = g.even_regions();
       int odd = g.odd_regions(); 
       bool mcurve = g.is_mcurve();
   }
```

