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

## Installation

#### Prerequisites
In order to compile the generators you need GCC, OpenMPI and [Google Sparsehash](https://github.com/sparsehash/sparsehash).
If you want to generate random delaunay graphs, you also need CGAL.
You can install these dependencies via your package manager:
```shell
sudo apt-get install gcc-11 g++-11 libopenmpi-dev libcgal-dev libsparsehash-dev 
sudo pacman -S gcc sparsehash openmpi cgal
```

#### Compiling 
To compile the code either run `compile.sh` or use the following instructions:
```shell
git submodule update --init --recursive
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

#### Output
By default our generators will output the generated graphs as a DIMACS edge list. 
Other output formats are available:

- `-f edge-list`: DIMACS edge list format (default)
- `-f binary-edge-list`: DIMACS binary edge list format
- `-f metis`: Metis undirected graph format (only works with undirected graph generators)
- `-f hmetis`: hMetis hypergraph format 

For a list of all output options, see `--help`.

#### As a Library 

To use KaGen as a library, simply add this repository as a Git submodule, 
include its CMake script using `add_subdirectory(extern/KaGen)` and link 
against the `KaGen::KaGen` target:

```cmake 
add_subdirectory(external/KaGen)
target_link_libraries(<target> PUBLIC KaGen::KaGen)
```

## Graph Models

### Erdos-Renyi Graphs G(n, m)
Generate a random graph using the Erdos-Renyi model G(n, m).
The graph can either be directed or undirected and can contain self-loops.

#### Parameters
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

#### Interface
```c++
KaGen gen(proc_rank, proc_size);
auto [edge_list_directed, vertex_range_directed] = gen.GenerateDirectedGNM(n, m, self_loops);
auto [edge_list_undirected, vertex_range_undirected] = gen.GenerateUndirectedGNM(n, m, self_loops);
```

---

### Erdos-Renyi Graphs G(n, p)
Generate a random graph using the Erdos-Renyi model G(n, p).
The graph can either be directed or undirected and can contain self-loops.

#### Parameters
```
mpirun -n <nproc> ./kagen <gnp_directed|gnp_undirected> 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -p <edge probability>
  [--self-loops]
  [-k <number of chunks>]
  [-S <seed>]
```

#### Interface
```c++
KaGen gen(proc_rank, proc_size);
auto [edge_list_directed, vertex_range_directed] = gen.GenerateDirectedGNP(n, p, self_loops);
auto [edge_list_undirected, vertex_range_undirected] = gen.GenerateUndirectedGNP(n, p, self_loops);
```

---

### Random Geometric Graphs RGG(n, r)
Generate a random graph using the random geometric graph model RGG(n, r).
NOTE: The number of processes must be a power of 2. 

#### Parameters
```
mpirun -n <nproc> ./kagen <rgg2d|rgg3d> 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -r <edge radius>
  [-k <number of chunks>]
  [-S <seed>]
```

#### Interface
```c++
KaGen gen(proc_rank, proc_size);
auto [edge_list_2d, vertex_range_2d] = gen.Generate2DRGG(n, r);
auto [edge_list_3d, vertex_range_3d] = gen.Generate3DRGG(n, r);
```

--- 

### Random Delaunay Graphs RDG(n)
Generate a random graph using the random Delaunay graph model RDG(n).
NOTE: The graph is generated with periodic boundary conditions to avoid long edges at the border.

#### Parameters
```
mpirun -n <nproc> ./kagen <rdg2d|rdg3d>
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  [-k <number of chunks>]
  [-S <seed>]
```

#### Interface
```c++
KaGen gen(proc_rank, proc_size);
auto [edge_list_2d, vertex_range_2d] = gen.Generate2DRDG(n);
auto [edge_list_3d, vertex_range_3d] = gen.Generate3DRDG(n);
```

---

### Barabassi-Albert Graphs BA(n, d)
Generate a random directed graph using the Barabassi-Albert graph model BA(n, d).
The graph may contain self-loops and multi edges.

#### Parameters
```
mpirun -n <nproc> ./kagen ba 
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -d <minimum degree for each vertex>
  [-k <number of chunks>]
  [-S <seed>]
```

#### Interface
```c++
KaGen gen(proc_rank, proc_size);
auto [edge_list, vertex_range] = gen.GenerateBA(n, md);
```

--- 

### Random Hyperbolic Graphs RHG(n, gamma, d)
Generate a two dimensional random graph using the random hyperbolic graph model RHG(n, gamma, d).

#### Parameters
```
mpirun -n <nproc> ./kagen rhg
  -n <number of vertices>
  [-N <number of vertices as a power of two>]
  -g <power-law exponent>
  -d <average vertex degree>
  [-k <number of chunks>]
  [-S <seed>]
```

#### Interface
```c++
KaGen gen(proc_rank, proc_size);
auto [edge_list, vertex_range] = gen.GenerateRHG(n, gamma, d);
```

--- 

**[License](/LICENSE):** 2-clause BSD
