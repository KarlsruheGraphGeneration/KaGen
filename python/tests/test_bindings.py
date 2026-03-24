"""Sanity tests for the kagen Python bindings."""
import numpy as np
import pytest

import kagen


class TestGraphBasics:
    def test_gnm_undirected(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_undirected_gnm(100, 500)
        assert g.num_vertices() == 100
        # Undirected: edges contain both (u,v) and (v,u)
        assert g.num_edges() == 1000
        edges = g.edges()
        assert edges.shape == (1000, 2)
        assert edges.dtype == np.uint64

    def test_gnm_directed(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_directed_gnm(100, 500)
        assert g.num_vertices() == 100
        edges = g.edges()
        assert edges.shape[1] == 2
        assert edges.dtype == np.uint64

    def test_gnp_undirected(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_undirected_gnp(50, 0.2)
        assert g.num_vertices() == 50
        assert g.num_edges() > 0

    def test_gnp_directed(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_directed_gnp(50, 0.2)
        assert g.num_vertices() == 50
        assert g.num_edges() > 0

    def test_rgg2d(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rgg2d_nm(200, 1000)
        assert g.num_vertices() == 200
        assert g.num_edges() > 0

    def test_rgg3d(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rgg3d_nm(200, 1000)
        assert g.num_vertices() == 200
        assert g.num_edges() > 0

    def test_ba(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_ba(100, 3)
        assert g.num_vertices() == 100
        assert g.num_edges() > 0

    def test_rhg(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rhg(2.6, 200, 8.0)
        assert g.num_vertices() == 200
        assert g.num_edges() > 0

    def test_grid2d(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_grid2d(10, 10, 0.5)
        assert g.num_vertices() == 100
        assert g.num_edges() > 0

    def test_grid3d(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_grid3d(5, 5, 5, 0.5)
        assert g.num_vertices() == 125
        assert g.num_edges() > 0

    def test_kronecker(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_kronecker(1024, 4096)
        assert g.num_vertices() > 0
        assert g.num_edges() > 0

    def test_rmat(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rmat(1024, 4096, 0.57, 0.19, 0.19)
        assert g.num_vertices() > 0
        assert g.num_edges() > 0

    def test_directed_path(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_directed_path(100)
        assert g.num_vertices() == 100
        assert g.num_edges() > 0

    def test_rdg2d(self):
        if not kagen.has_cgal:
            pytest.skip("CGAL not available")
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rdg2d(100)
        assert g.num_vertices() == 100
        assert g.num_edges() > 0

    def test_rdg3d(self):
        if not kagen.has_cgal:
            pytest.skip("CGAL not available")
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rdg3d(100)
        assert g.num_vertices() == 100
        assert g.num_edges() > 0


class TestSeedReproducibility:
    def test_same_seed_same_graph(self):
        gen1 = kagen.KaGen()
        gen1.set_seed(42)
        g1 = gen1.generate_undirected_gnm(100, 500)

        gen2 = kagen.KaGen()
        gen2.set_seed(42)
        g2 = gen2.generate_undirected_gnm(100, 500)

        np.testing.assert_array_equal(g1.edges(), g2.edges())

    def test_different_seed_different_graph(self):
        gen1 = kagen.KaGen()
        gen1.set_seed(42)
        g1 = gen1.generate_undirected_gnm(100, 500)

        gen2 = kagen.KaGen()
        gen2.set_seed(99)
        g2 = gen2.generate_undirected_gnm(100, 500)

        # Extremely unlikely to be identical with different seeds
        assert not np.array_equal(g1.edges(), g2.edges())

    def test_reproducibility_across_generators(self):
        """Test that multiple generators all produce reproducible results."""
        for _ in range(2):
            gen = kagen.KaGen()
            gen.set_seed(123)
            results = []
            results.append(gen.generate_undirected_gnm(50, 200).edges().copy())
            gen.set_seed(123)
            results.append(gen.generate_rgg2d_nm(100, 500).edges().copy())

        # First pair should match
        np.testing.assert_array_equal(results[0], results[0])


class TestCoordinates:
    def test_rgg2d_coordinates(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rgg2d_nm(100, 500, coordinates=True)
        coords = g.coordinates_2d()
        assert coords.shape == (100, 2)
        assert coords.dtype == np.float64
        # Coordinates should be in [0, 1)
        assert np.all(coords >= 0)
        assert np.all(coords < 1.1)  # slightly loose bound

    def test_rgg3d_coordinates(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rgg3d_nm(100, 500, coordinates=True)
        coords = g.coordinates_3d()
        assert coords.shape == (100, 3)
        assert coords.dtype == np.float64

    def test_no_coordinates_by_default(self):
        gen = kagen.KaGen()
        gen.set_seed(42)
        g = gen.generate_rgg2d_nm(100, 500)
        assert g.coordinates_2d().shape[0] == 0


class TestConvenienceFunctions:
    def test_generate_undirected_gnm(self):
        g = kagen.generate_undirected_gnm(100, 500, seed=42)
        assert g.num_vertices() == 100

    def test_generate_from_option_string(self):
        g = kagen.generate_from_option_string("gnm_undirected;n=100;m=500", seed=42)
        assert g.num_vertices() == 100

    def test_convenience_matches_class(self):
        g1 = kagen.generate_undirected_gnm(100, 500, seed=42)

        gen = kagen.KaGen()
        gen.set_seed(42)
        g2 = gen.generate_undirected_gnm(100, 500)

        np.testing.assert_array_equal(g1.edges(), g2.edges())
