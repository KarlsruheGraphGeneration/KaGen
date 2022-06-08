# Communication-free Graph Generators 

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

Additionally, if you use the Barabassi-Albert generator, we ask that you cite the paper:

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
* The random delaunay graph generators also require CGAL. 

You can install these dependencies via your package manager:

```shell
# Ubuntu, Debian 
apt-get install gcc-11 g++-11 libopenmpi-dev libcgal-dev libsparsehash-dev 

# Arch Linux, Manjaro
pacman -S gcc sparsehash openmpi cgal

# Fedora 
dnf install gcc openmpi sparsehash-devel CGAL-devel

# macOS using Homebrew 
brew install gcc open-mpi google-sparsehash cgal

# macOS using MacPorts 
port install gcc11 openmpi sparsehash cgal5
```

## Building KaGen 
To compile the code either run `compile.sh` or use the following instructions:

```shell
git submodule update --init --recursive
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make kagen kagen_library -j
```

## Running KaGen 

After building KaGen, the standalone application is located at `build/app/kagen`. 
A list of all command line options is available using the `./kagen --help` option. 
To view the options of a specific graph generator, use:

```shell 
./kagen <gnm_undirected|gnm_directed|gnp_undirected|gnp_directed|rgg2d|rgg3d|grid2d|grid3d|rdg2d|rdg3d|rhg|ba|kronecker> --help
```

By default, the generated graph is written to a single file `out.edgelist` (`-o` option) in DIMACS edge list format (`-f` option).
Other output formats are available:

- `-f edge-list`: DIMACS edge list format (default)
- `-f binary-edge-list`: DIMACS binary edge list format
- `-f metis`: Metis graph format
- `-f hmetis`: hMetis hypergraph format 
- `-f dot`: GraphViz dot file (add `-C` to include vertex coordinates for 2D graph generators)
- `-f coordinates`: Text file containing vertex coordinates 

If you want each PE to write its edges to a seperate file, use the `--distributed-output` flag.

## Using the KaGen Library

The KaGen library is located at `build/library/libkagen_library.a` (use `-DBUILD_SHARED_LIBS=On` to build a shared library instead).
Alternatively, you can use KaGen by adding this repository as a Git submodule to your project and including it in your CMake configuration:

```cmake 
add_subdirectory(external/KaGen)
target_link_libraries(<your-target> PUBLIC KaGen::KaGen)
```

## Graph Models

### Erdos-Renyi Graphs G(n, m)
Generate a random graph using the Erdos-Renyi model G(n, m).
The graph can either be directed or undirected and can contain self-loops.

#### Application
```
mpirun -n <nproc> ./kagen <gnm_directed|gnm_undirected> 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -m <number of edges>
  [-M <number of edges as a power of two>]
  [--self-loops]
  [-k <number of chunks>]
  [-S <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);
auto [edge_list_directed, vertex_range_directed] = gen.GenerateDirectedGNM(n, m, self_loops);
auto [edge_list_undirected, vertex_range_undirected] = gen.GenerateUndirectedGNM(n, m, self_loops);
```

---

### Erdos-Renyi Graphs G(n, p)
Generate a random graph using the Erdos-Renyi model G(n, p).
The graph can either be directed or undirected and can contain self-loops.

#### Application
```
mpirun -n <nproc> ./kagen <gnp_directed|gnp_undirected> 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -p <edge probability>
  [--self-loops]
  [-k <number of chunks>]
  [-S <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);
auto [edge_list_directed, vertex_range_directed] = gen.GenerateDirectedGNP(n, p, self_loops);
auto [edge_list_undirected, vertex_range_undirected] = gen.GenerateUndirectedGNP(n, p, self_loops);
```

---

### Random Geometric Graphs RGG(n, r)
Generate an undirected random graph using the random geometric graph model RGG(n, r).

**Note:** This generator requires the number of PEs to be a power of 2.

**Note:** This generator is parameterized by the number of nodes in the graph and its edge radius. 
Either parameter can be omitted in favor of the desired number of edges, in which case the omitted 
parameter is approximated such that the expected number of edges matches the desired number of edges.

#### Application
```
mpirun -n <nproc> ./kagen <rgg2d|rgg3d> 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -r <edge radius>
  -m <number of edges>                     # only if -n or -r are omitted
  [-M <number of edges as a power of two>] # only if -n or -r are omitted
  [-k <number of chunks>] 
  [-S <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);
auto [edge_list_2d, vertex_range_2d] = gen.GenerateRGG2D(n, r);
auto [edge_list_2d, vertex_range_2d] = gen.GenerateRGG2DNM(n, m); // deduce r s.t. E[# edges] = m
auto [edge_list_2d, vertex_range_2d] = gen.GenerateRGG2DMR(m, r); // deduce n s.t. E[# edges] = m
auto [edge_list_3d, vertex_range_3d] = gen.GenerateRGG3D(n, r);
auto [edge_list_3d, vertex_range_3d] = gen.GenerateRGG3DNM(n, m); // deduce r s.t. E[# edges] = m
auto [edge_list_3d, vertex_range_3d] = gen.GenerateRGG3DMR(m, r); // deduce n s.t. E[# edges] = m
```

--- 

### Random Delaunay Graphs RDG(n)
Generate an undirected random graph using the random Delaunay graph model RDG(n).

**Note:** The graph is generated with periodic boundary conditions to avoid long edges at the border.

#### Application
```
mpirun -n <nproc> ./kagen <rdg2d|rdg3d>
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  [-k <number of chunks>]
  [-S <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);
auto [edge_list_2d, vertex_range_2d] = gen.GenerateRDG2D(n);
auto [edge_list_3d, vertex_range_3d] = gen.GenerateRDG3D(n);
```

---

### Random Grid Graphs GRID(x, y [, z])
Generate an undirected random grid graph. 

#### Application 
```
mpirun -n <nproc> ./kagen <grid2d|grid3d>
  -x <width of grid>
  [-X <width of grid as a power of two>]
  -y <height of grid>
  [-Y <height of grid as a power of two>]
  -z <depth of grid (grid3d only)>
  [-Z <depth of grid as a power of two (grid3d only)>]
  [-k <number of cunks>]
  [-S <seed>]
```

#### Library 
```c++
KaGen gen(MPI_COMM_WORLD);
auto [edge_list_2d, vertex_range_2d] = gen.GenerateGrid2D(x, y);
auto [edge_list_3d, vertex_range_3d] = gen.GenerateGrid3D(x, y, z);
```

---

### Barabassi-Albert Graphs BA(n, d)
Generate a random directed graph using the Barabassi-Albert graph model BA(n, d).
The graph may contain self-loops and multi edges.

#### Application
```
mpirun -n <nproc> ./kagen ba 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -d <minimum degree for each vertex>
  [-k <number of chunks>]
  [-S <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);
auto [edge_list, vertex_range] = gen.GenerateBA(n, md);
```

--- 

### Random Hyperbolic Graphs RHG(n, gamma, d)
Generate a two dimensional undirected random graph using the random hyperbolic graph model RHG(n, gamma, d).

#### Application
```
mpirun -n <nproc> ./kagen rhg
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -g <power-law exponent>
  -d <average vertex degree>
  [-k <number of chunks>]
  [-S <seed>]
```

#### Library
```c++
KaGen gen(MPI_COMM_WORLD);
auto [edge_list, vertex_range] = gen.GenerateRHG(n, gamma, d);
```

--- 

**[License](/LICENSE):** 2-clause BS
