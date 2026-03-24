#!/usr/bin/env python3
"""
KaGen Python Example: Generate a graph from every available model.

Usage:
    python python_example.py
"""
import kagen
import numpy as np

SEED = 42

print(f"KaGen Python v{kagen.__version__}")
print(f"CGAL support: {kagen.has_cgal}")
print("=" * 78)


def show(name, graph, extra=""):
    n = graph.num_vertices()
    m = graph.num_edges()
    edges = graph.edges()
    if m > 0:
        max_degree = np.bincount(edges[:, 0].astype(np.int64)).max()
    else:
        max_degree = 0
    line = f"  {name:<35s}  n={n:>6,}  m={m:>8,}  max_deg={max_degree:>4}"
    if extra:
        line += f"  {extra}"
    print(line)


# ─── Erdos-Renyi Models ────────────────────────────────────────────────────
print("\nErdos-Renyi Models:")

gen = kagen.KaGen()
gen.set_seed(SEED)

g = gen.generate_undirected_gnm(1000, 5000)
show("G(n,m) undirected", g)

g = gen.generate_directed_gnm(1000, 5000)
show("G(n,m) directed", g)

g = gen.generate_undirected_gnp(500, 0.02)
show("G(n,p) undirected", g, f"p=0.02")

g = gen.generate_directed_gnp(500, 0.02)
show("G(n,p) directed", g, f"p=0.02")


# ─── Random Geometric Graphs ──────────────────────────────────────────────
print("\nRandom Geometric Graphs:")

g = gen.generate_rgg2d_nm(2000, 10000, coordinates=True)
coords = g.coordinates_2d()
show("RGG 2D (n,m)", g, f"coords={coords.shape}")

g = gen.generate_rgg2d(1000, 0.05, coordinates=True)
show("RGG 2D (n,r)", g, f"r=0.05")

g = gen.generate_rgg3d_nm(2000, 10000, coordinates=True)
coords = g.coordinates_3d()
show("RGG 3D (n,m)", g, f"coords={coords.shape}")

g = gen.generate_rgg3d(1000, 0.1, coordinates=True)
show("RGG 3D (n,r)", g, f"r=0.1")


# ─── Random Delaunay Graphs (requires CGAL) ───────────────────────────────
if kagen.has_cgal:
    print("\nRandom Delaunay Graphs:")

    g = gen.generate_rdg2d(1000, periodic=False, coordinates=True)
    show("RDG 2D", g, f"coords={g.coordinates_2d().shape}")

    g = gen.generate_rdg2d(1000, periodic=True, coordinates=True)
    show("RDG 2D (periodic)", g)

    g = gen.generate_rdg3d(500, coordinates=True)
    show("RDG 3D", g, f"coords={g.coordinates_3d().shape}")
else:
    print("\nRandom Delaunay Graphs: SKIPPED (no CGAL)")


# ─── Random Hyperbolic Graphs ─────────────────────────────────────────────
print("\nRandom Hyperbolic Graphs:")

g = gen.generate_rhg(2.6, 5000, 10.0, coordinates=True)
show("RHG (gamma,n,d)", g, f"gamma=2.6, d=10.0")

g = gen.generate_rhg_nm(2.6, 5000, 25000)
show("RHG (gamma,n,m)", g, f"gamma=2.6")


# ─── Barabasi-Albert ──────────────────────────────────────────────────────
print("\nBarabasi-Albert (Preferential Attachment):")

g = gen.generate_ba(5000, 3)
show("BA (n,d)", g, f"min_degree=3")

g = gen.generate_ba(5000, 3, directed=True)
show("BA directed", g, f"min_degree=3")


# ─── Grid Graphs ──────────────────────────────────────────────────────────
print("\nGrid Graphs:")

g = gen.generate_grid2d(50, 50, 0.5, coordinates=True)
show("Grid 2D (50x50, p=0.5)", g)

g = gen.generate_grid2d(50, 50, 1.0, periodic=True)
show("Grid 2D (periodic, p=1)", g)

g = gen.generate_grid3d(10, 10, 10, 0.5)
show("Grid 3D (10x10x10, p=0.5)", g)

g = gen.generate_grid3d(10, 10, 10, 1.0, periodic=True)
show("Grid 3D (periodic, p=1)", g)


# ─── Kronecker / R-MAT ───────────────────────────────────────────────────
print("\nKronecker / R-MAT:")

g = gen.generate_kronecker(1024, 8192)
show("Kronecker", g)

g = gen.generate_rmat(1024, 8192, 0.57, 0.19, 0.19)
show("R-MAT (a=.57,b=.19,c=.19)", g)

g = gen.generate_rmat(1024, 8192, 0.57, 0.19, 0.19, directed=True)
show("R-MAT directed", g)


# ─── Path Graph ───────────────────────────────────────────────────────────
print("\nPath / Cycle:")

g = gen.generate_directed_path(1000)
show("Directed path", g)

g = gen.generate_directed_path(1000, periodic=True)
show("Directed cycle (periodic)", g)


# ─── Option String Interface ──────────────────────────────────────────────
print("\nOption String Interface:")

g = gen.generate_from_option_string("gnm_undirected;n=1000;m=5000")
show("gnm_undirected;n=1000;m=5000", g)

g = gen.generate_from_option_string("rgg2d;n=2000;m=10000")
show("rgg2d;n=2000;m=10000", g)

g = gen.generate_from_option_string("rhg;n=5000;gamma=2.6;avg_degree=10")
show("rhg;n=5000;gamma=2.6;d=10", g)


# ─── Convenience Functions ────────────────────────────────────────────────
print("\nConvenience Functions (module-level):")

g = kagen.generate_undirected_gnm(1000, 5000, seed=42)
show("kagen.generate_undirected_gnm()", g)

g = kagen.generate_rgg2d_nm(2000, 10000, seed=42, coordinates=True)
show("kagen.generate_rgg2d_nm()", g, f"coords={g.coordinates_2d().shape}")

g = kagen.generate_rhg(2.6, 5000, 10.0, seed=42)
show("kagen.generate_rhg()", g)


# ─── CSR Representation ──────────────────────────────────────────────────
print("\nCSR Representation:")

gen_csr = kagen.KaGen()
gen_csr.set_seed(SEED)
gen_csr.use_csr_representation()
g = gen_csr.generate_undirected_gnm(1000, 5000)
xadj = g.xadj()
adjncy = g.adjncy()
print(f"  GNM CSR: xadj.shape={xadj.shape}, adjncy.shape={adjncy.shape}")
print(f"  Vertex 0 has {xadj[1]-xadj[0]} neighbors: {adjncy[xadj[0]:xadj[1]][:5]}...")


# ─── Edge & Vertex Weights ───────────────────────────────────────────────
print("\nEdge & Vertex Weights:")

gen_w = kagen.KaGen()
gen_w.set_seed(SEED)
gen_w.configure_edge_weight_generation("hashing_based", 1, 100)
gen_w.configure_vertex_weight_generation("uniform_random", 1, 50)
g = gen_w.generate_undirected_gnm(500, 2000)
ew = g.edge_weights()
vw = g.vertex_weights()
print(f"  Edge weights:   count={len(ew)}, range=[{ew.min()}, {ew.max()}]")
print(f"  Vertex weights: count={len(vw)}, range=[{vw.min()}, {vw.max()}]")


print("\n" + "=" * 78)
print("All models generated successfully!")
