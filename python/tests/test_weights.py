"""Tests for edge and vertex weight generation."""
import numpy as np
import pytest

import kagen


class TestEdgeWeights:
    def test_hashing_based(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        gen.configure_edge_weight_generation("hashing_based", 1, 100)
        g = gen.generate_undirected_gnm(100, 500)
        ew = g.edge_weights()
        assert len(ew) == g.num_edges()
        assert ew.dtype == np.int64
        assert np.all(ew >= 1)
        assert np.all(ew < 100)

    def test_uniform_random(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        gen.configure_edge_weight_generation("uniform_random", 10, 50)
        g = gen.generate_undirected_gnm(100, 500)
        ew = g.edge_weights()
        assert len(ew) == g.num_edges()
        assert np.all(ew >= 10)
        assert np.all(ew < 50)

    def test_no_weights_by_default(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_undirected_gnm(100, 500)
        assert len(g.edge_weights()) == 0

    def test_weight_reproducibility(self):
        for _ in range(2):
            gen = kagen.KaGen()
            gen.set_seed(42)
            gen.configure_edge_weight_generation("hashing_based", 1, 1000)
            g = gen.generate_undirected_gnm(100, 500)
            ew = g.edge_weights()
        # The last iteration should produce the same weights
        gen2 = kagen.KaGen()
        gen2.set_seed(42)
        gen2.configure_edge_weight_generation("hashing_based", 1, 1000)
        g2 = gen2.generate_undirected_gnm(100, 500)
        np.testing.assert_array_equal(ew, g2.edge_weights())


class TestVertexWeights:
    def test_uniform_random(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        gen.configure_vertex_weight_generation("uniform_random", 1, 100)
        g = gen.generate_undirected_gnm(100, 500)
        vw = g.vertex_weights()
        assert len(vw) == g.num_vertices()
        assert vw.dtype == np.int64
        assert np.all(vw >= 1)
        assert np.all(vw < 100)

    def test_no_weights_by_default(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_undirected_gnm(100, 500)
        assert len(g.vertex_weights()) == 0
