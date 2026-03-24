"""
KaGen: Communication-free Massively Distributed Graph Generators.

Python bindings for the KaGen library providing fast, deterministic graph
generation for various random graph models.
"""

from ._kagen import KaGen, Graph, has_cgal

__version__ = "1.2.1"
__all__ = [
    "KaGen",
    "Graph",
    "has_cgal",
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
