"""Tests for GenerateFromOptionString."""
import numpy as np

import kagen
from conftest import sorted_edges


def test_option_string_gnm():
    gen1 = kagen.KaGen()
    gen1.set_seed(42)
    g1 = gen1.generate_undirected_gnm(1024, 4096)

    gen2 = kagen.KaGen()
    gen2.set_seed(42)
    g2 = gen2.generate_from_option_string("gnm_undirected;n=1024;m=4096")

    np.testing.assert_array_equal(sorted_edges(g1), sorted_edges(g2))


def test_option_string_rgg2d():
    gen1 = kagen.KaGen()
    gen1.set_seed(42)
    g1 = gen1.generate_rgg2d_nm(512, 2048)

    gen2 = kagen.KaGen()
    gen2.set_seed(42)
    g2 = gen2.generate_from_option_string("rgg2d;n=512;m=2048")

    np.testing.assert_array_equal(sorted_edges(g1), sorted_edges(g2))


def test_option_string_rhg():
    gen1 = kagen.KaGen()
    gen1.set_seed(42)
    g1 = gen1.generate_rhg(2.6, 1024, 8.0)

    gen2 = kagen.KaGen()
    gen2.set_seed(42)
    g2 = gen2.generate_from_option_string("rhg;n=1024;gamma=2.6;avg_degree=8.0")

    np.testing.assert_array_equal(sorted_edges(g1), sorted_edges(g2))


def test_option_string_ba():
    gen1 = kagen.KaGen()
    gen1.set_seed(42)
    g1 = gen1.generate_ba(512, 3)

    gen2 = kagen.KaGen()
    gen2.set_seed(42)
    g2 = gen2.generate_from_option_string("ba;n=512;min_degree=3")

    np.testing.assert_array_equal(sorted_edges(g1), sorted_edges(g2))
