# Communication-free Graph Generators (+ others)

This is the code to accompany our eponymous paper: *Funke, D., Lamm, S., Sanders, P., Schulz, C., Strash, D. and von Looz, M., 2017. Communication-free Massively Distributed Graph Generation. arXiv preprint arXiv:1710.07565.*
You can find a freely accessible online version [in the arXiv](https://arxiv.org/abs/1710.07565).

If you use this library in the context of an academic publication, we ask that you cite our paper:

```bibtex
@inproceedings{funke2017communication,
  title={Communication-free Massively Distributed Graph Generation},
  author={Funke, Daniel and Lamm, Sebastian and Sanders, Peter and Schulz, Christian and Strash, Darren and von Looz, Moritz},
  booktitle={2018 {IEEE} International Parallel and Distributed Processing Symposium, {IPDPS} 2018, Vancouver, BC, Canada, May 21 -- May 25, 2018},
  year={2018},
}
```

Additionally, if you use the Barabassi-Albert generator, we ask that you cite the [paper](https://arxiv.org/abs/1602.07106):

```bibtex
@article{sanders2016generators,
  title={Scalable generation of scale-free graphs},
  journal={Information Processing Letters},
  volume={116},
  number={7},
  pages={489 -- 491},
  year={2016},
  author={Sanders, Peter and Schulz, Christian},
}
```

If you use the R-MAT generator, we ask that you cite the [paper](https://www.cambridge.org/core/journals/network-science/article/linear-work-generation-of-rmat-graphs/68A0DDA58A7B84E9B3ACA2DBB123A16C):

```bibtex
@article{HubSan2020RMAT, 
  title={Linear Work Generation of {R-MAT} Graphs}, 
  volume={8}, 
  number={4}, 
  journal={Network Science}, 
  publisher={Cambridge University Press}, 
  author={H{\"u}bschle-Schneider, Lorenz and Sanders, Peter}, 
  year={2020}, 
  pages={543 -- 550},
}
```

## Introduction 

Network generators serve as a tool to alleviate the need for synthethic instances with controllable parameters by algorithm developers and researchers. 
However, many generators fail to provide instances on a massive scale due to their sequential nature or resource constraints.

In our work, we present novel generators for a variety of network models commonly found in practice.
By making use of pseudorandomization and divide-and-conquer schemes, our generators follow a communication-free paradigm.
The resulting generators are often embarrassingly parallel and have a near optimal scaling behavior.
This allows us to generate instances of up to 2^43 vertices and 2^47 edges in less than 22 minutes on 32768 cores.
Therefore, our generators allow new graph families to be used on an unprecedented scale.

## Requirements 

In order to compile the generators, you require: 

* A modern, C++17-ready compiler such as `g++` version 9 or higher or `clang` version 11 or higher. 
  * Note: Apple Clang is **not** supported. 
* OpenMPI
* [Google Sparsehash](https://github.com/sparsehash/sparsehash)
* CGAL (optional, only required for the Delaunay generators)

You can install these dependencies via your package manager:

```shell
# Ubuntu, Debian 
apt-get install gcc-12 g++-12 libopenmpi-dev libcgal-dev libsparsehash-dev 

# Arch Linux, Manjaro
pacman -S gcc sparsehash openmpi cgal

# Fedora 
dnf install gcc openmpi sparsehash-devel CGAL-devel

# macOS using Homebrew 
brew install gcc open-mpi google-sparsehash cgal

# macOS using MacPorts 
port install gcc12 openmpi sparsehash cgal5
```

## Building KaGen 

To compile the code either run `compile.sh` or use the following instructions:

```shell
git submodule update --init --recursive
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

## Running KaGen 

After building KaGen, the standalone application is located at `build/app/KaGen`. 
A list of all command line options is available using the `./KaGen --help` option. 
To view the options of a specific graph generator, use:

```shell 
./KaGen <gnm-undirected|gnm-directed|gnp-undirected|gnp-directed|rgg2d|rgg3d|grid2d|grid3d|rdg2d|rdg3d|rhg|ba|kronecker|rmat> --help
```

By default, the generated graph is written to a single file `out` (`-o` option) in DIMACS edge list format (`-f` option).
Other output formats include:

- `-f edgelist`: DIMACS edge list format (default)
- `-f binary-edgelist`: DIMACS binary edge list format, use `--32` to write the file with 32 bit data types 
- `-f metis`: Metis graph format
- `-f hmetis`: hMetis hypergraph format; **note:** KaGen still generates a graph, i.e., every hyperedge will contain two pins
- `-f dot`: GraphViz dot file (add `-C` to include vertex coordinates for 2D graph generators)
- `-f coordinates`: Text file containing vertex coordinates 
- `-f parhip`: Binary graph format used by [ParHIP](https://github.com/KaHIP/KaHIP)
- `-f xtrapulp`: Binary graph format used by [XtraPuLP](https://github.com/HPCGraphAnalysis/PuLP), use `--32` to write the file with 32 bit data types

Experimental output formats include:

- `-f experimental/hmetis-ep`: hMetis hypergraph format, but the graph is transformed s.t. a partition of the hypergraph is an edge partition of the generated graph
- `-f experimental/freight-netl`: hypergraph format used by FREIGHT; **note:** KaGen still generates a graph, i.e., every hyperedge will contain two pins
- `-f experimental/freight-netl-ep`: hypergraph format used by FREIGHT, but the graph is transformed s.t. a partition of the hypergraph is an edge partition of the generated graph

One graph can be stored in multiple formats by passing the `-f` repeatedly, e.g., `-o out -f metis -f coordinates` will write two files `out.metis` and `out.xyz`.
If you want each PE to write its edges to a seperate file, use the `--distributed-output` flag.

## Using the KaGen Library

The KaGen library is located at `build/library/libkagen.a` (use `-DBUILD_SHARED_LIBS=On` to build a shared library instead) and can be used in C++ and C projects.
If you are using CMake, you can use KaGen by adding this repository as a Git submodule to your project and including it in your CMake configuration:

```cmake 
add_subdirectory(external/KaGen)
target_link_libraries(<your-target> PUBLIC KaGen::KaGen)
```

Alternatively, you can use `FetchContent`: 

```cmake 
include(FetchContent)
FetchContent_Declare(KaGen 
  GIT_REPOSITORY https://github.com/sebalamm/KaGen.git 
  GIT_TAG master)
FetchContent_MakeAvailable(KaGen)
set_property(DIRECTORY "${KaGen_SOURCE_DIR}" PROPERTY EXCLUDE_FROM_ALL YES) # optional

target_link_libraries(<your-target> PUBLIC KaGen::KaGen)
```

Examples on how to use the C and C++ interfaces are available in the `examples/` directory.
The examples given below only show the C++ interface.

**Note**: Instead of calling the library functions listed below, you can also use `KaGen::GenerateFromOptionString()` 
to pass the generator options as a string (documentation is available in `kagen/kagen.h`).

The library functions return the generated graph as an instance of type `kagen::Graph`. 
By default, the graph is represented as an edge list, i.e., a vector `kagen::Graph::edges[]` containing pairs of vertices.
To generate a graph in compressed sparse row (CSR) format, call `kagen::KaGen::UseCSRRepresentation()` before generating the graph. 
Then, access the graph via `kagen::Graph::xadj[]` and `kagen::Graph::adjncy[]`.

## General Graph Format

Unless noted otherwise, KaGen generates **simple**, **undirected** graphs, i.e.,
graphs without self-loops, without multi-edges and where for every edge (u, v),
there is also a reverse edge (v, u).

When using KaGen in a distributed setting, each PE owns an equally sized range of consecutive vertices.
An edge is owned by the PE that owns its tail vertex.
Thus, an edge (u, v) is in the edge list of the PE that owns vertex u, while the reverse edge 
(v, u) is in the edge list of the PE owning vertex v.

## Communication-free Graph Generators

### Erdos-Renyi Graphs with Fixed Number of Edges
Generate a random Erdos-Renyi graph with a fixed number of edges.
The graph can either be directed or undirected and can contain self-loops.

#### Application
```
mpirun -n <nproc> ./KaGen <gnm-directed|gnm-undirected> 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -m <number of edges>
  [-M <number of edges as a power of two>]
  [--self-loops]
  [-k <number of chunks>]
  [-s <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph_directed = gen.GenerateDirectedGNM(n, m, self_loops = false);
Graph graph_undirected = gen.GenerateUndirectedGNM(n, m, self_loops = false);
```

---

### Erdos-Renyi graphs with Fixed Edge Probability
Generate a random Erdos-Renyi graph with a fixed edge probability.
The graph can either be directed or undirected and can contain self-loops.

#### Application
```
mpirun -n <nproc> ./KaGen <gnp_directed|gnp_undirected> 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -p <edge probability>
  [--self-loops]
  [-k <number of chunks>]
  [-s <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph_directed = gen.GenerateDirectedGNP(n, p, self_loops = false);
Graph graph_undirected = gen.GenerateUndirectedGNP(n, p, self_loops = false);
```

---

### Random Geometric Graphs
Generate an undirected random geometric graph.

**Note:** This generator is parameterized by the number of vertices in the graph and its edge radius. 
Either parameter can be omitted in favor of the desired number of edges, in which case the omitted 
parameter is approximated such that the expected number of edges matches the desired number of edges.

#### Application
```
mpirun -n <nproc> ./KaGen <rgg2d|rgg3d> 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -r <edge radius>
  -m <number of edges>                     # only if -n or -r are omitted
  [-M <number of edges as a power of two>] # only if -n or -r are omitted
  [-k <number of chunks>] 
  [-s <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateRGG2D(n, r, coordinates = false);
Graph graph = gen.GenerateRGG2D_NM(n, m, coordinates = false); // deduce r s.t. E[# edges] = m
Graph graph = gen.GenerateRGG2D_MR(m, r, coordinates = false); // deduce n s.t. E[# edges] = m

Graph graph = gen.GenerateRGG3D(n, r, coordinates = false);
Graph graph = gen.GenerateRGG3D_NM(n, m, coordinates = false); // deduce r s.t. E[# edges] = m
Graph graph = gen.GenerateRGG3D_MR(m, r, coordinates = false); // deduce n s.t. E[# edges] = m
```

--- 

### Random Delaunay Graphs
Generate an undirected random delaunay graph.

**Note:** The graph can be generated with periodic boundary conditions to avoid long edges at the border using the `-p` flag. 
However, this can yield unexpected results when using less than 9 PEs (2D) / 27 PEs (3D) to generate the graph.

#### Application
```
mpirun -n <nproc> ./KaGen <rdg2d|rdg3d>
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  [--periodic]
  [-s <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateRDG2D(n, periodic, coordinates = false);
Graph graph = gen.GenerateRDG2D_M(m, periodic, coordinates = false);

Graph graph = gen.GenerateRDG3D(n, coordinates = false);
Graph graph = gen.GenerateRDG3D_M(m, coordinates = false);
```

---

### Random Grid Graphs
Generate an undirected random grid graph. 

#### Application 
```
mpirun -n <nproc> ./KaGen <grid2d|grid3d>
  -x <width of grid>
  [-X <width of grid as a power of two>]
  -y <height of grid>
  [-Y <height of grid as a power of two>]
  -z <depth of grid (grid3d only)>
  [-Z <depth of grid as a power of two (grid3d only)>]
  -p <edge probability>
  -m <number of edges>                     # only if -p is omitted
  [-M <number of edges as a power of two>] # only if -p is omitted
  [--periodic]
  [-k <number of cunks>]
  [-s <seed>]
```

#### Library 
```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateGrid2D(x, y, p, periodic, coordinates = false);
Graph graph = gen.GenerateGrid2D_N(n, p, periodic, coordinates = false); // x, y = sqrt(n)
Graph graph = gen.GenerateGrid2D_NM(n, m, periodic, coordinates = false); // x, y = sqrt(n)

Graph graph = gen.GenerateGrid3D(x, y, z, p, periodic, coordinates = false);
Graph graph = gen.GenerateGrid3D_N(n, p, periodic, coordinates = false); // x, y, z = cbrt(n) 
Graph graph = gen.GenerateGrid3D_NM(n, m, periodic, coordinates = false); // x, y, z = cbrt(n) 
```

--- 

### Random Hyperbolic Graphs 
Generate a two dimensional undirected random hyperbolic graph.

**Note:** On x86 systems, the generator can use 64 bit or 80 bit floating point numbers.
This can be controlled explicitly by using the `--hp-floats` or `--no-hp-floats` flags. 
If neither flag is set, KaGen switches to 80 bit precision automatically if the generated graph has more than 2^29 vertices.

**Note:** Due to floating point inaccuracies, this generator performs communication in a post-processing step.

#### Application
```
mpirun -n <nproc> ./KaGen rhg
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -g <power-law exponent>
  -d <average vertex degree>
  [-k <number of chunks>]
  [--hp-floats]
  [--no-hp-floats]
  [-s <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateRHG(gamma, n, d, coordinates = false);
Graph graph = gen.GenerateRHG_NM(gamma, n, m, coordinates = false); // deduce d s.t. E[# edges] = m
Graph graph = gen.GenerateRHG_MD(gamma, m, d, coordinates = false); // deduce n s.t. E[# edges] = m
```

## Non-communication-free Graph Generators 

Since the original publication, several other graph generators have been integrated into the KaGen framework. 

### Barabassi-Albert Graphs 

Generate a random Barabassi-Albert graph.
The graph may contain self-loops and multi edges.

#### Application
```
mpirun -n <nproc> ./KaGen ba 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -d <minimum degree for each vertex>
  [--directed]
  [--self-loops]
  [-k <number of chunks>]
  [-s <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateBA(n, d, directed = false, self_loops = false);
Graph graph = gen.GenerateBA_NM(n, m, directed = false, self_loops = false);
Graph graph = gen.GenerateBA_MD(m, d, directed = false, self_loops = false);
```

---

### R-MAT Graphs
Generate a random RMAT graph.

Each PE generates a random R-MAT graph with n vertices and m/\<nproc\> edges.
Afterwards, the vertices are assigned to PEs round-robin style and edges are distributed accordingly.

#### Application
```
mpirun -n <nproc> ./KaGen rmat 
  -n <number of vertices> # should be a power of two
  [-N <number of vertices as a power of two>]
  -m <number of edges>
  [-M <number of edges as a power of two>]
  -a <probability for an edge to land in block a>
  -b <probability for an edge to land in block b>
  -c <probability for an edge to land in block c>
  [--directed]
  [--self-loops]
  [-s <seed>]
```

#### Library 
```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateRMAT(n, m, a, b, c, directed = false, self_loops = false);
```

---

### Kronecker Graphs 
Generate a random Kronecker graph.

Each PE generates a random Kronecker graph with n vertices and m/\<nproc\> edges.
Afterwards, the vertices are assigned to PEs round-robin style and edges are distributed accordingly.

#### Application 
```
mpirun -n <nproc> ./KaGen kronecker 
  -n <number of vertices> # should be a power of two 
  [-N <number of vertices as a power of  two>]
  -m <number of edges> 
  [-M <number of edges as a power of two>]
  [--directed]
  [--self-loops]
  [-s <seed>]
```

#### Library 

```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateKronecker(n, m, directed = false, self_loops = false);
```

## Static Graph Generators

### Image Graph Generator
Generates a graph based on an input image.
Each pixel is represented by a vertex with edges to its neighboring vertices.
The image has to be converted to KARGB format first (a simple binary file containing the uncompressed R, G, B channels of the image) by 
using the `img2kargb` or `upsb2kargb` tool shipped with KaGen.

#### Application
```
mpirun -n <nproc> ./KaGen image
  --filename=<path to kargb file>
  [--weight-model=<l2, inv-l2, inv-ratio>]
  [--weight-multiplier=1]
  [--weight-offset=0]
  [--min-weight-threshold=1]
  [--max-weight-threshold=inf]
  [--neighborhood=<4, 8, 24>]
  [--max-grid-x=<...>]
  [--max-grid-y=<...>]
  [--grid-x=<...>]
  [--grid-y=<...>]
  [--cols-per-pe=<...>]
  [--rows-per-pe=<...>]
```

#### Library 

```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateFromOptionString("image;filename=<...>;...");
```

--- 

### File Graph Generator
Pseudo-generator that loads a static graph from disk.
Can be used to convert input formats to output format, or to load static graphs when using KaGen as a library.

#### Application 
```
mpirun -n <nproc> ./KaGen file
  --filename=<path to graph>
  --input-format=<metis|parhip>
  [--distribution=<balance-vertices|balance-edges>]
```

#### Library 

```c++
KaGen gen(MPI_COMM_WORLD);

Graph graph = gen.GenerateFromOptionString("file;filename=<...>;input_format=<...>;distribution=<...>");
```

## Tools

Tools can be installed via `cmake --install build --component tools`. The following tools are included: 

```shell
# graphstats: compute some basic statistics for the given graphs
mpirun ./app/tools/graphstats <path to graph(s), ...>
  [-f <format, e.g., metis, parhip, plain-edgelist>]

# chkgraph: validate a graph file in any supported input format
mpirun -n <nproc> ./app/tools/chkgraph <path to graph>
  [-f <format, e.g., metis, parhip, plain-edgelist>] 
  [--64bits]                  # allow 64 bit weights and IDs
  [--self-loops]              # allow self loops
  [--directed]                # allow directed graphs (i.e., not all reverse edges are present)
  [--multi-edges]             # allow multi edges
  [--negative-edge-weights]   # allow negative edge weights
  [--negative-vertex-weights] # allow negative vertex weights
  
# pangraph: convert a graph file between supported formats in external memory
./app/tools/pangraph --input-format=<...> --input-filename=<...> --output-format=<...> --output-filename=<...>
  [-C <num chunks = 1>]       # split the graph into <num chunks> chunks; only one chunk has to fit into internal memory at a time
  [-T <tmp directory = /tmp>] # directory to be used for temporary files (requires free space roughly the size of the input graph)
  [--remove-self-loops]       # remove any self-loops during convertion
  [--add-reverse-edges]       # make all edges undirected by adding potentially missing reverse edges
  [--sort-edges]              # sort the outgoing edges by destination vertex ID
  [-n <num vertices>]         # provide the number of vertices in the graph -- currently only used for the plain-edgelist input format
```

---

**[License](/LICENSE):** 2-clause BS
