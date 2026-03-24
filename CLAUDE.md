# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

KaGen (Karlsruhe Graph Generation) is a communication-free massively distributed graph generator library. It generates synthetic graphs from various random graph models (Erdos-Renyi, Random Geometric, Hyperbolic, Barabasi-Albert, Kronecker, R-MAT, Grid, etc.) at massive scale. Written in C++ with Python bindings via pybind11.

## Build Commands

```bash
# Full build (init submodules + cmake + compile)
./compile.sh

# Manual build
git submodule update --init --recursive
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

# Debug build
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel

# Build with tests enabled
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DKAGEN_BUILD_TESTS=ON
cmake --build build --parallel

# Build with Python bindings (no MPI, for local Python use)
cmake -B build -DKAGEN_BUILD_PYTHON=ON -DKAGEN_NOMPI=ON
cmake --build build --parallel
```

Key CMake options: `KAGEN_NOMPI` (pseudo-MPI for single-process), `KAGEN_USE_CGAL` (Delaunay generators), `KAGEN_BUILD_APPS`, `KAGEN_BUILD_TOOLS`, `KAGEN_BUILD_EXAMPLES`, `KAGEN_WARNINGS_ARE_ERRORS`.

## Testing

```bash
# C++ tests (GoogleTest via KaTestrophe, MPI-aware)
cd build && ctest --output-on-failure

# Python tests
pytest python/tests

# Run a single C++ test
cd build && ctest -R <test_name> --output-on-failure
```

CI runs both GCC and Clang in Debug and Release modes with `-DKAGEN_WARNINGS_ARE_ERRORS=ON`.

## Code Formatting

Uses clang-format (config in `.clang-format`). Column limit: 120. C++20 standard (set in `.clangd`).

## Architecture

- **`kagen/`** — Core library
  - **`kagen.h` / `kagen.cpp`** — Public C++ and C API. The `KaGen` class is the main entry point; `Graph` holds results (edge list or CSR, with optional weights/coordinates).
  - **`generators/`** — Each subdirectory implements a graph model (gnm, gnp, geometric, hyperbolic, barabassi, kronecker, rmat, grid, path, image, file). All inherit from `Generator` base class with `GenerateEdgeList()` and `GenerateCSR()` methods.
  - **`io/`** — I/O format handlers (METIS, DIMACS, DOT, ParHIP, etc.)
  - **`edgeweight_generators/`** and **`vertexweight_generators/`** — Post-generation weight assignment
  - **`context.h`** — Configuration holder passed through the system
  - **`sampling/`** — Sampling utilities used by generators
  - **`tools/`** — Internal utilities (RNG, geometry, hashing)
- **`python/`** — pybind11 bindings (`src/_kagen.cpp`) and Python module (`kagen/__init__.py` with convenience functions)
- **`app/`** — CLI application (`KaGen.cpp`) and tools (graphstats, chkgraph, pangraph)
- **`tests/`** — GoogleTest-based tests using KaTestrophe (MPI-aware test framework)
- **`extlib/`** — Bundled: libmorton (Z-order curves), xxHash
- **`nompi/`** — Pseudo-MPI stub for `KAGEN_NOMPI` builds

## Key Types (kagen/kagen.h)

- `SInt` = `unsigned long long` (vertex/edge indices)
- `SSInt` = `long long` (weights)
- `HPFloat` = `long double`, `LPFloat` = `double`
- `Edgelist` = `vector<pair<SInt, SInt>>`
- Graph representations: `EDGE_LIST` or `CSR` (compressed sparse row)

## Dependencies

Required: MPI (OpenMPI), C++17+ compiler (GCC 9+, Clang 11+; Apple Clang not supported). Optional: CGAL (Delaunay), Google Sparsehash, xxHash, Intel MKL. Python bindings need pybind11 and numpy.
