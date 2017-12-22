# Communication-free Graph Generators

This is the code to accompany our eponymous paper: *Funke, D., Lamm, S., Sanders, P., Schulz, C., Strash, D. and von Looz, M., 2017. Communication-free Massively Distributed Graph Generation. arXiv preprint arXiv:1710.07565.*
You can find a freely accessible online version [in the arXiv](https://arxiv.org/abs/1710.07565).

If you use this library in the context of an academic publication, we ask that you cite our paper:
```bibtex
@article{funke2017communication,
  title={Communication-free Massively Distributed Graph Generation},
  author={Funke, Daniel and Lamm, Sebastian and Sanders, Peter and Schulz, Christian and Strash, Darren and von Looz, Moritz},
  journal={arXiv preprint arXiv:1710.07565},
  year={2017}
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

## Installation

#### Prerequisites
To compile the generators you need [Google Sparsehash](https://github.com/sparsehash/sparsehash).

Compile instructions:
```
  mkdir build
  cd build
  cmake ..
  make
```

By default our generators will output the generated graphs in the DIMACS format.
For performance measurements this can be turned of by disabling the `-DOUTPUT_EDGES` flag.

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

#### Example
Generate a directed G(n,m) graph with 2^20 vertices and 2^22 edges with self-loops on 16 processors and write it to tmp
```
mpirun -n 16 ./build/kagen -gen gnm_directed -n 20 -m 22 -self_loops -output tmp
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

#### Example
Generate a directed G(n,p) graph with 2^20 vertices and an edge probability of 0.001 with self-loops on 16 processors and write it to tmp
```
mpirun -n 16 ./build/kagen -gen gnp_directed -n 20 -p 0.001 -self_loops -output tmp
```

---

#### Random Geometric Graphs RGG(n,r)
Generate a random graph using the random geometric graph model RGG(n,r).
Graphs will always be undirected and can be either two- or three-dimensional.
NOTE: Use a square (cubic) number of chunks/processes for the two-dimensional (three-dimensional) generator.
##### Parameters
```
-gen <rgg_2d|rgg_3d>
-n <number of vertices as a power of two>
-r <radius for vertices to be connected> (r <= 1.0)
-k <number of chunks>
-seed <seed for PRNGs>
-output <output file>
```

##### Example
Generate a three dimensional RGG(n,r) graph with 2^20 vertices and a radius of 0.00275 on 16 processors and write it to tmp
```
mpirun -n 16 ./build/kagen -gen rgg_3d -n 20 -r 0.00275 -output tmp
```

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

##### Example
Generate a three dimensional RDG(n,r) graph with 2^20 vertices on 16 processors and write it to tmp
```
mpirun -n 16 ./build/kagen -gen rgg_2d -n 20 -output tmp
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

##### Example
Generate a BA(n,d) graph with 2^20 vertices and a minimum degree of 4 on 16 processors and write it to tmp
```
mpirun -n 16 ./build/kagen -gen ba -n 20 -d 4 -output tmp
```

--- 

#### Random Hyperbolic Graphs RHG(n,d)
Generate a two dimensional random graph using the random hyperbolic graph model RHG(n,gamma,d)$
##### Parameters
```
-gen rhg
-n <number of vertices as a power of two>
-d <average degree> 
-gamma <power-law exponent> 
-k <number of chunks>
-seed <seed for PRNGs>
-output <output file>
```

##### Example
Generate a two dimensional RHG(n,r) graph with 2^20 vertices and an average degree of 8 with a power-law exponent of 2.2 on 16 processors and write it to tmp
```
mpirun -n 16 ./build/kagen -gen rhg -n 20 -d 8 -gamma 2.2 -output tmp
```

--- 

**[License](/LICENSE):** 2-clause BSD
