"""
Tests that verify Python bindings produce exactly the same graphs as
the KaGen CLI (mpirun -np 1).
"""
import numpy as np
import pytest

import kagen
from conftest import run_kagen_cli, sorted_edges


SEED = 42


# Each entry: (test_id, cli_generator, cli_args, python_method, python_kwargs)
GENERATOR_CONFIGS = [
    (
        "gnm_undirected",
        "gnm-undirected",
        ["-n", "1024", "-m", "4096"],
        "generate_undirected_gnm",
        {"n": 1024, "m": 4096},
    ),
    (
        "gnm_directed",
        "gnm-directed",
        ["-n", "1024", "-m", "4096"],
        "generate_directed_gnm",
        {"n": 1024, "m": 4096},
    ),
    (
        "gnp_undirected",
        "gnp-undirected",
        ["-n", "256", "-p", "0.1"],
        "generate_undirected_gnp",
        {"n": 256, "p": 0.1},
    ),
    (
        "gnp_directed",
        "gnp-directed",
        ["-n", "256", "-p", "0.1"],
        "generate_directed_gnp",
        {"n": 256, "p": 0.1},
    ),
    (
        "rgg2d",
        "rgg2d",
        ["-n", "512", "-m", "2048"],
        "generate_rgg2d_nm",
        {"n": 512, "m": 2048},
    ),
    (
        "rgg3d",
        "rgg3d",
        ["-n", "512", "-m", "2048"],
        "generate_rgg3d_nm",
        {"n": 512, "m": 2048},
    ),
    (
        "ba",
        "ba",
        ["-n", "512", "-d", "3"],
        "generate_ba",
        {"n": 512, "d": 3},
    ),
    (
        "rhg",
        "rhg",
        ["-g", "2.6", "-n", "1024", "--avg-deg", "8.0"],
        "generate_rhg",
        {"gamma": 2.6, "n": 1024, "d": 8.0},
    ),
    (
        "grid2d",
        "grid2d",
        ["-x", "32", "-y", "32", "-p", "0.5"],
        "generate_grid2d",
        {"grid_x": 32, "grid_y": 32, "p": 0.5},
    ),
    (
        "grid3d",
        "grid3d",
        ["-x", "8", "-y", "8", "-z", "8", "-p", "0.5"],
        "generate_grid3d",
        {"grid_x": 8, "grid_y": 8, "grid_z": 8, "p": 0.5},
    ),
    (
        "kronecker",
        "kronecker",
        ["-n", "1024", "-m", "4096"],
        "generate_kronecker",
        {"n": 1024, "m": 4096},
    ),
    (
        "rmat",
        "rmat",
        ["-n", "1024", "-m", "4096", "-a", "0.57", "-b", "0.19", "-c", "0.19"],
        "generate_rmat",
        {"n": 1024, "m": 4096, "a": 0.57, "b": 0.19, "c": 0.19},
    ),
    (
        "path_directed",
        "path-directed",
        ["-n", "256"],
        "generate_directed_path",
        {"n": 256},
    ),
]

# CGAL-dependent generators
CGAL_GENERATOR_CONFIGS = [
    (
        "rdg2d",
        "rdg2d",
        ["-n", "256"],
        "generate_rdg2d",
        {"n": 256, "periodic": False},
    ),
    (
        "rdg3d",
        "rdg3d",
        ["-n", "256"],
        "generate_rdg3d",
        {"n": 256},
    ),
]


@pytest.mark.cli
@pytest.mark.parametrize(
    "test_id,cli_gen,cli_args,py_method,py_kwargs",
    GENERATOR_CONFIGS,
    ids=[c[0] for c in GENERATOR_CONFIGS],
)
def test_cli_reproducibility(kagen_binary, test_id, cli_gen, cli_args, py_method, py_kwargs):
    """Python must produce the exact same sorted edge list as the CLI."""
    # Generate via CLI
    cli_edges = run_kagen_cli(kagen_binary, cli_gen, cli_args, seed=SEED)

    # Generate via Python
    gen = kagen.KaGen()
    gen.set_seed(SEED)
    method = getattr(gen, py_method)
    graph = method(**py_kwargs)
    py_edges = sorted_edges(graph)

    np.testing.assert_array_equal(
        cli_edges,
        py_edges,
        err_msg=f"Edge mismatch for {test_id}: CLI has {len(cli_edges)} edges, Python has {len(py_edges)}",
    )


@pytest.mark.cli
@pytest.mark.parametrize(
    "test_id,cli_gen,cli_args,py_method,py_kwargs",
    CGAL_GENERATOR_CONFIGS,
    ids=[c[0] for c in CGAL_GENERATOR_CONFIGS],
)
def test_cli_reproducibility_cgal(kagen_binary, test_id, cli_gen, cli_args, py_method, py_kwargs):
    """CGAL generators: Python must produce the exact same sorted edge list as the CLI."""
    if not kagen.has_cgal:
        pytest.skip("CGAL not available")

    cli_edges = run_kagen_cli(kagen_binary, cli_gen, cli_args, seed=SEED)

    gen = kagen.KaGen()
    gen.set_seed(SEED)
    method = getattr(gen, py_method)
    graph = method(**py_kwargs)
    py_edges = sorted_edges(graph)

    np.testing.assert_array_equal(
        cli_edges,
        py_edges,
        err_msg=f"Edge mismatch for {test_id}: CLI has {len(cli_edges)} edges, Python has {len(py_edges)}",
    )
