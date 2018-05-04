# Communication-free Graph Generators [![Build Status](https://travis-ci.com/sebalamm/KaGen_Release.svg?token=xPJBtRVg2pEbu9qYXkx1&branch=master)](https://travis-ci.com/sebalamm/KaGen_Release) 

This is the code to accompany our eponymous paper: *Funke, D., Lamm, S., Sanders, P., Schulz, C., Strash, D. and von Looz, M., 2017. Communication-free Massively Distributed Graph Generation. arXiv preprint arXiv:1710.07565.*
You can find a freely accessible online version [in the arXiv](https://arxiv.org/abs/1710.07565).

If you use this library in the context of an academic publication, we ask that you cite our paper:
```bibtex
@inproceedings{funke2017communication,
  title={Communication-free Massively Distributed Graph Generation},
  author={Funke, Daniel and Lamm, Sebastian and Sanders, Peter and Schulz, Christian and Strash, Darren and von Looz, Moritz},
  booktitle={2018 {IEEE} International Parallel and Distributed Processing Symposium, {IPDPS} 2018, Vancouver, BC, Canada, May 21 - May 25, 2018},
  year={2018},
  intype={to appear in}
}
```

Additionally, if you use the Barabassi-Albert generator, we ask that you cite the paper:
```bibtex
@article{sanders2016generators,
  title = "Scalable generation of scale-free graphs",
  journal = "Information Processing Letters",
  volume = "116",
  number = "7",
  pages = "489 - 491",
  year = "2016",
  author = "Peter Sanders and Christian Schulz",
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
In order to compile the generators you need g++-7, OpenMPI, CGAL and [Google Sparsehash](https://github.com/sparsehash/sparsehash).
If you haven't installed these dependencies, please do so via your package manager.
```
  sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
  sudo apt-get -qq update
  sudo apt-get install gcc-7 g++-7 libopenmpi-dev libcgal-dev libcgal-qt5-dev libsparsehash-dev 
```

#### Compiling 
To compile the code either run `compile.sh` or use the following instruction
```
  git submodule update --init --recursive
  mkdir build
  cd build
  cmake ..
  make
```

By default our generators will output the generated graphs in the DIMACS format.
For performance measurements this can be turned of by disabling the `-DOUTPUT_EDGES` flag.
Additional flags for varying the output can be found in `CMakeLists.txt`.

## Graph Models

### Erdos-Renyi Graphs G(n,m)
Generate a random graph using the Erdos-Renyi model G(n,m).
The graph can either be directed or undirected and can contain self-loops.

#### Parameters
```
-gen <gnm_directed|gnm_undirected>
-n <number of vertices as a power of two>
-m <number of edges as a power of two>
-k <number of chunks> 
-seed <seed for PRNGs>
-output <output file>
-self_loops 
```

#### Interface
```
KaGen gen(proc_rank, proc_size);
auto edge_list_directed = gen.GenerateDirectedGNM(n, m, k, seed, output, self_loops);
auto edge_list_undirected = gen.GenerateUndirectedGNM(n, m, k, seed, output, self_loops);
```

#### Command Line Example
Generate a directed G(n,m) graph with 2^20 vertices and 2^22 edges with self-loops on 16 processors and write it to tmp
```
mpirun -n 16 ./build/app/kagen -gen gnm_directed -n 20 -m 22 -self_loops -output tmp
```

---

### Erdos-Renyi Graphs G(n,p)
Generate a random graph using the Erdos-Renyi model G(n,p).
The graph can either be directed or undirected and can contain self-loops.

#### Parameters
```
-gen <gnp_directed|gnp_undirected>
-n <number of vertices as a power of two>
-p <edge probability>
-k <number of chunks> 
-seed <seed for PRNGs>
-output <output file>
-self_loops 
```

#### Interface
```
KaGen gen(proc_rank, proc_size);
auto edge_list_directed = gen.GenerateDirectedGNP(n, p, k, seed, output, self_loops);
auto edge_list_undirected = gen.GenerateUndirectedGNP(n, p, k, seed, output, self_loops);
```

#### Command Line Example
Generate a directed G(n,p) graph with 2^20 vertices and an edge probability of 0.001 with self-loops on 16 processors and write it to tmp
```
mpirun -n 16 ./build/app/kagen -gen gnp_directed -n 20 -p 0.001 -self_loops -output tmp
```
kagen
--- 

#### Random Delaunay Graphs RDG(n)
Generate a random graph using the random Delaunay graph model RDG(n).
NOTE: Use a square (cubic) number of chunks/processes for the two-dimensional (three-dimensional) generator.
##### Parameters
```
-gen <rdg_2d|rdg_3d>
-n <number of vertices as a power of two>
-k <number of chunks>
-seed <seed for PRNGs>
-output <output file>
```

#### Interface
```
KaGen gen(proc_rank, proc_size);
auto edge_list_2d = gen.Generate2DRDG(n, k, seed, output);
auto edge_list_3d = gen.Generate3DRGG(n, k, seed, output);
```

##### Command Line Example
Generate a three dimensional RDG(n,r) graph with 2^20 vertices on 16 processors and write it to tmp
```
mpirun -n 16 ./build/app/kagen -gen rdg_3d -n 20 -output tmp
```

--- 

#### Barabassi-Albert Graphs BA(n,d)
Generate a random graph using the Barabassi-Albert graph model BA(n,d)
##### Parameters
```
-gen ba
-n <number of vertices as a power of two>
-md <min degree for each vertex> 
-k <number of chunks>
-seed <seed for PRNGs>
-output <output file>
```

#### Interface
```
KaGen gen(proc_rank, proc_size);
auto edge_list = gen.GenerateBA(n, md, seed, output);
```

##### Command Line Example
Generate a BA(n,d) graph with 2^20 vertices and a minimum degree of 4 on 16 processors and write it to tmp
```
mpirun -n 16 ./build/app/kagen -gen ba -n 20 -md 4 -output tmp
```

--- 

#### Random Hyperbolic Graphs RHG(n,gamma,d)
Generate a two dimensional random graph using the random hyperbolic graph model RHG(n,gamma,d)
##### Parameters
```
-gen rhg
-n <number of vertices as a power of two>
-gamma <power-law exponent> 
-d <average degree> 
-k <number of chunks>
-seed <seed for PRNGs>
-output <output file>
```

#### Interface
```
KaGen gen(proc_rank, proc_size);
auto edge_list = gen.GenerateRHG(n, gamma, d, seed, output);
```

##### Command Line Example
Generate a two dimensional RHG(n,r) graph with 2^20 vertices and an average degree of 8 with a power-law exponent of 2.2 on 16 processors and write it to tmp
```
mpirun -n 16 ./build/app/ -gen rhg -n 20 -d 8 -gamma 2.2 -output tmp
```

--- 

**[License](/LICENSE):** kagen2-clause BSD
