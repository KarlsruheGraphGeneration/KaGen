"""
KaGen: Communication-free Massively Distributed Graph Generators.

Python bindings for the KaGen library providing fast, deterministic graph
generation for various random graph models.
"""

from ._kagen import KaGen, Graph, has_cgal

__version__ = "1.2.9"


def describe():
    """Return a comprehensive description of the kagen module for agentic AI systems.

    This function provides all information needed to programmatically use kagen:
    available generators, their parameters, return types, and example calls.

    Returns:
        str: A structured description of the module.
    """
    return '''# kagen - Communication-free Graph Generators

## Overview
KaGen generates synthetic graphs from various random graph models.
All generators are deterministic given a seed. Returns numpy arrays.
CGAL support (for Delaunay generators): ''' + str(has_cgal) + '''

## Graph Object
All generators return a `kagen.Graph` with these methods:
  - graph.edges()           -> np.ndarray, shape (num_edges, 2), dtype=uint64
  - graph.xadj()            -> np.ndarray, CSR row offsets, dtype=uint64
  - graph.adjncy()          -> np.ndarray, CSR column indices, dtype=uint64
  - graph.vertex_weights()  -> np.ndarray, dtype=int64 (empty if not configured)
  - graph.edge_weights()    -> np.ndarray, dtype=int64 (empty if not configured)
  - graph.coordinates_2d()  -> np.ndarray, shape (n, 2), dtype=float64
  - graph.coordinates_3d()  -> np.ndarray, shape (n, 3), dtype=float64
  - graph.vertex_range()    -> tuple (first_vertex, first_invalid_vertex)
  - graph.num_vertices()    -> int
  - graph.num_edges()       -> int
  - graph.sort_edgelist()   -> None (sorts edges in-place)

Note: undirected graphs store both (u,v) and (v,u), so num_edges() == 2 * logical_edges.

## Convenience Functions (module-level)
All accept `seed=<int>` for reproducibility. Quick single-call interface.

### Erdos-Renyi (random edges)
  kagen.generate_undirected_gnm(n, m, seed=None, self_loops=False)
  kagen.generate_directed_gnm(n, m, seed=None, self_loops=False)
  kagen.generate_undirected_gnp(n, p, seed=None, self_loops=False)
  kagen.generate_directed_gnp(n, p, seed=None, self_loops=False)
  # n=num_vertices, m=num_edges, p=edge_probability

### Random Geometric (edges between nearby points)
  kagen.generate_rgg2d(n, r, seed=None, coordinates=False)
  kagen.generate_rgg2d_nm(n, m, seed=None, coordinates=False)
  kagen.generate_rgg3d(n, r, seed=None, coordinates=False)
  kagen.generate_rgg3d_nm(n, m, seed=None, coordinates=False)
  # r=radius threshold, coordinates=True to get vertex positions

### Random Delaunay (triangulation of random points, requires CGAL)
  kagen.generate_rdg2d(n, seed=None, periodic=False, coordinates=False)
  kagen.generate_rdg3d(n, seed=None, coordinates=False)

### Random Hyperbolic (power-law degree distribution)
  kagen.generate_rhg(gamma, n, d, seed=None, coordinates=False)
  kagen.generate_rhg_nm(gamma, n, m, seed=None, coordinates=False)
  # gamma=power-law exponent (typically 2.1-3.0), d=average_degree

### Barabasi-Albert (preferential attachment, scale-free)
  kagen.generate_ba(n, d, seed=None, directed=False, self_loops=False)
  kagen.generate_ba_nm(n, m, seed=None, directed=False, self_loops=False)
  # d=minimum_degree for each new vertex

### Grid (regular lattice with random edge removal)
  kagen.generate_grid2d(grid_x, grid_y, p, seed=None, periodic=False, coordinates=False)
  kagen.generate_grid3d(grid_x, grid_y, grid_z, p, seed=None, periodic=False, coordinates=False)
  # p=probability each grid edge is included (1.0 = full grid)

### Kronecker / R-MAT (recursive matrix models)
  kagen.generate_kronecker(n, m, seed=None, directed=False, self_loops=False)
  kagen.generate_rmat(n, m, a, b, c, seed=None, directed=False, self_loops=False)
  # a,b,c=block probabilities (a+b+c < 1, d=1-a-b-c). Typical: a=0.57, b=c=0.19

### Path / Cycle
  kagen.generate_directed_path(n, seed=None, permute=False, periodic=False)
  # periodic=True makes it a cycle

### Option String (alternative interface)
  kagen.generate_from_option_string("gnm_undirected;n=1000;m=5000", seed=42)
  # Syntax: "model;key=value;key=value"

## Class-based Interface (for multiple generations or advanced config)
  gen = kagen.KaGen()
  gen.set_seed(42)
  gen.use_csr_representation()          # default is edge list
  gen.configure_edge_weight_generation("hashing_based", 1, 100)
  gen.configure_vertex_weight_generation("uniform_random", 1, 50)
  graph = gen.generate_undirected_gnm(1000, 5000)
  # Weight generator types: "hashing_based", "uniform_random", "euclidean_distance"

## Example Calls

```python
import kagen
import numpy as np

# Simple Erdos-Renyi graph
g = kagen.generate_undirected_gnm(n=1000, m=5000, seed=42)
edges = g.edges()  # shape (10000, 2) -- undirected stores both directions

# Random geometric graph with coordinates
g = kagen.generate_rgg2d_nm(n=2000, m=10000, seed=1, coordinates=True)
coords = g.coordinates_2d()  # shape (2000, 2)

# Power-law graph (social network-like)
g = kagen.generate_rhg(gamma=2.6, n=10000, d=10.0, seed=42)

# Scale-free graph (preferential attachment)
g = kagen.generate_ba(n=5000, d=3, seed=42)

# Grid graph
g = kagen.generate_grid2d(grid_x=100, grid_y=100, p=1.0, seed=42)

# CSR format with weights
gen = kagen.KaGen()
gen.set_seed(42)
gen.use_csr_representation()
gen.configure_edge_weight_generation("uniform_random", 1, 100)
g = gen.generate_undirected_gnm(500, 2000)
xadj, adjncy, weights = g.xadj(), g.adjncy(), g.edge_weights()
```
'''


__all__ = [
    "KaGen",
    "Graph",
    "has_cgal",
    "describe",
    # Convenience functions
    "generate_from_option_string",
    "generate_directed_gnm",
    "generate_undirected_gnm",
    "generate_directed_gnp",
    "generate_undirected_gnp",
    "generate_rgg2d",
    "generate_rgg2d_nm",
    "generate_rgg2d_mr",
    "generate_rgg3d",
    "generate_rgg3d_nm",
    "generate_rgg3d_mr",
    "generate_rdg2d",
    "generate_rdg2d_m",
    "generate_rdg3d",
    "generate_rdg3d_m",
    "generate_ba",
    "generate_ba_nm",
    "generate_ba_md",
    "generate_rhg",
    "generate_rhg_nm",
    "generate_rhg_md",
    "generate_grid2d",
    "generate_grid2d_n",
    "generate_grid2d_nm",
    "generate_grid3d",
    "generate_grid3d_n",
    "generate_grid3d_nm",
    "generate_directed_path",
    "generate_kronecker",
    "generate_rmat",
]


def _make_gen(seed=None):
    gen = KaGen()
    if seed is not None:
        gen.set_seed(seed)
    return gen


def generate_from_option_string(options, seed=None):
    """Generate a graph from an option string, e.g. 'gnm_undirected;n=100;m=500'."""
    return _make_gen(seed).generate_from_option_string(options)


def generate_directed_gnm(n, m, seed=None, self_loops=False):
    """Generate a directed Erdos-Renyi G(n,m) graph."""
    return _make_gen(seed).generate_directed_gnm(n, m, self_loops)


def generate_undirected_gnm(n, m, seed=None, self_loops=False):
    """Generate an undirected Erdos-Renyi G(n,m) graph."""
    return _make_gen(seed).generate_undirected_gnm(n, m, self_loops)


def generate_directed_gnp(n, p, seed=None, self_loops=False):
    """Generate a directed Erdos-Renyi G(n,p) graph."""
    return _make_gen(seed).generate_directed_gnp(n, p, self_loops)


def generate_undirected_gnp(n, p, seed=None, self_loops=False):
    """Generate an undirected Erdos-Renyi G(n,p) graph."""
    return _make_gen(seed).generate_undirected_gnp(n, p, self_loops)


def generate_rgg2d(n, r, seed=None, coordinates=False):
    """Generate a 2D random geometric graph."""
    return _make_gen(seed).generate_rgg2d(n, r, coordinates)


def generate_rgg2d_nm(n, m, seed=None, coordinates=False):
    """Generate a 2D random geometric graph with n vertices and m edges."""
    return _make_gen(seed).generate_rgg2d_nm(n, m, coordinates)


def generate_rgg2d_mr(m, r, seed=None, coordinates=False):
    """Generate a 2D random geometric graph with m edges and radius r."""
    return _make_gen(seed).generate_rgg2d_mr(m, r, coordinates)


def generate_rgg3d(n, r, seed=None, coordinates=False):
    """Generate a 3D random geometric graph."""
    return _make_gen(seed).generate_rgg3d(n, r, coordinates)


def generate_rgg3d_nm(n, m, seed=None, coordinates=False):
    """Generate a 3D random geometric graph with n vertices and m edges."""
    return _make_gen(seed).generate_rgg3d_nm(n, m, coordinates)


def generate_rgg3d_mr(m, r, seed=None, coordinates=False):
    """Generate a 3D random geometric graph with m edges and radius r."""
    return _make_gen(seed).generate_rgg3d_mr(m, r, coordinates)


def generate_rdg2d(n, seed=None, periodic=False, coordinates=False):
    """Generate a 2D random Delaunay graph (requires CGAL)."""
    return _make_gen(seed).generate_rdg2d(n, periodic, coordinates)


def generate_rdg2d_m(m, seed=None, periodic=False, coordinates=False):
    """Generate a 2D random Delaunay graph with m edges (requires CGAL)."""
    return _make_gen(seed).generate_rdg2d_m(m, periodic, coordinates)


def generate_rdg3d(n, seed=None, coordinates=False):
    """Generate a 3D random Delaunay graph (requires CGAL)."""
    return _make_gen(seed).generate_rdg3d(n, coordinates)


def generate_rdg3d_m(m, seed=None, coordinates=False):
    """Generate a 3D random Delaunay graph with m edges (requires CGAL)."""
    return _make_gen(seed).generate_rdg3d_m(m, coordinates)


def generate_ba(n, d, seed=None, directed=False, self_loops=False):
    """Generate a Barabasi-Albert preferential attachment graph."""
    return _make_gen(seed).generate_ba(n, d, directed, self_loops)


def generate_ba_nm(n, m, seed=None, directed=False, self_loops=False):
    """Generate a Barabasi-Albert graph with n vertices and m edges."""
    return _make_gen(seed).generate_ba_nm(n, m, directed, self_loops)


def generate_ba_md(m, d, seed=None, directed=False, self_loops=False):
    """Generate a Barabasi-Albert graph with m edges and minimum degree d."""
    return _make_gen(seed).generate_ba_md(m, d, directed, self_loops)


def generate_rhg(gamma, n, d, seed=None, coordinates=False):
    """Generate a random hyperbolic graph."""
    return _make_gen(seed).generate_rhg(gamma, n, d, coordinates)


def generate_rhg_nm(gamma, n, m, seed=None, coordinates=False):
    """Generate a random hyperbolic graph with n vertices and m edges."""
    return _make_gen(seed).generate_rhg_nm(gamma, n, m, coordinates)


def generate_rhg_md(gamma, m, d, seed=None, coordinates=False):
    """Generate a random hyperbolic graph with m edges and average degree d."""
    return _make_gen(seed).generate_rhg_md(gamma, m, d, coordinates)


def generate_grid2d(grid_x, grid_y, p, seed=None, periodic=False, coordinates=False):
    """Generate a 2D grid graph with edge probability p."""
    return _make_gen(seed).generate_grid2d(grid_x, grid_y, p, periodic, coordinates)


def generate_grid2d_n(n, p, seed=None, periodic=False, coordinates=False):
    """Generate a 2D grid graph with n vertices and edge probability p."""
    return _make_gen(seed).generate_grid2d_n(n, p, periodic, coordinates)


def generate_grid2d_nm(n, m, seed=None, periodic=False, coordinates=False):
    """Generate a 2D grid graph with n vertices and m edges."""
    return _make_gen(seed).generate_grid2d_nm(n, m, periodic, coordinates)


def generate_grid3d(grid_x, grid_y, grid_z, p, seed=None, periodic=False, coordinates=False):
    """Generate a 3D grid graph with edge probability p."""
    return _make_gen(seed).generate_grid3d(grid_x, grid_y, grid_z, p, periodic, coordinates)


def generate_grid3d_n(n, p, seed=None, periodic=False, coordinates=False):
    """Generate a 3D grid graph with n vertices and edge probability p."""
    return _make_gen(seed).generate_grid3d_n(n, p, periodic, coordinates)


def generate_grid3d_nm(n, m, seed=None, periodic=False, coordinates=False):
    """Generate a 3D grid graph with n vertices and m edges."""
    return _make_gen(seed).generate_grid3d_nm(n, m, periodic, coordinates)


def generate_directed_path(n, seed=None, permute=False, periodic=False):
    """Generate a directed path graph."""
    return _make_gen(seed).generate_directed_path(n, permute, periodic)


def generate_kronecker(n, m, seed=None, directed=False, self_loops=False):
    """Generate a Kronecker graph."""
    return _make_gen(seed).generate_kronecker(n, m, directed, self_loops)


def generate_rmat(n, m, a, b, c, seed=None, directed=False, self_loops=False):
    """Generate an R-MAT graph."""
    return _make_gen(seed).generate_rmat(n, m, a, b, c, directed, self_loops)
