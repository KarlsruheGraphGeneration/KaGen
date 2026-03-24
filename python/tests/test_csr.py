"""Tests that CSR representation is consistent with edge list."""
import numpy as np
import pytest

import kagen


def edges_from_csr(xadj, adjncy, vertex_range):
    """Reconstruct edge list from CSR representation."""
    first_vertex = vertex_range[0]
    edges = []
    for i in range(len(xadj) - 1):
        u = first_vertex + i
        for j in range(xadj[i], xadj[i + 1]):
            v = adjncy[j]
            edges.append((u, v))
    if not edges:
        return np.empty((0, 2), dtype=np.uint64)
    return np.array(edges, dtype=np.uint64)


@pytest.mark.parametrize("generator,kwargs", [
    ("generate_undirected_gnm", {"n": 100, "m": 500}),
    ("generate_directed_gnm", {"n": 100, "m": 500}),
    ("generate_rgg2d_nm", {"n": 200, "m": 1000}),
    ("generate_ba", {"n": 100, "d": 3}),
    ("generate_rhg", {"gamma": 2.6, "n": 200, "d": 8.0}),
    ("generate_grid2d", {"grid_x": 10, "grid_y": 10, "p": 0.5}),
])
def test_csr_matches_edgelist(generator, kwargs):
    """CSR and edge list representations must produce the same edges."""
    seed = 42

    # Edge list
    gen_el = kagen.KaGen()
    gen_el.set_seed(seed)
    gen_el.use_edge_list_representation()
    g_el = getattr(gen_el, generator)(**kwargs)
    el_edges = g_el.edges()

    # CSR
    gen_csr = kagen.KaGen()
    gen_csr.set_seed(seed)
    gen_csr.use_csr_representation()
    g_csr = getattr(gen_csr, generator)(**kwargs)
    csr_edges = edges_from_csr(g_csr.xadj(), g_csr.adjncy(), g_csr.vertex_range())

    # Sort both and compare
    if el_edges.shape[0] > 0:
        el_sorted = el_edges[np.lexsort((el_edges[:, 1], el_edges[:, 0]))]
    else:
        el_sorted = el_edges

    if csr_edges.shape[0] > 0:
        csr_sorted = csr_edges[np.lexsort((csr_edges[:, 1], csr_edges[:, 0]))]
    else:
        csr_sorted = csr_edges

    np.testing.assert_array_equal(el_sorted, csr_sorted)
