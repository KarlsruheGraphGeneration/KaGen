import os
import subprocess
import tempfile

import numpy as np
import pytest


def find_kagen_binary():
    """Find the KaGen CLI binary."""
    candidates = [
        os.path.join(os.path.dirname(__file__), "..", "..", "build-mpi", "app", "KaGen"),
        os.path.join(os.path.dirname(__file__), "..", "..", "build", "app", "KaGen"),
    ]
    for path in candidates:
        path = os.path.abspath(path)
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    return None


@pytest.fixture(scope="session")
def kagen_binary():
    """Path to the KaGen CLI binary built with MPI."""
    binary = find_kagen_binary()
    if binary is None:
        pytest.skip("KaGen CLI binary not found (build with MPI first)")
    return binary


def run_kagen_cli(binary, generator, args, seed=42):
    """Run KaGen CLI and return sorted edge list as numpy array.

    Returns edges as 0-indexed numpy array of shape (M, 2).
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        outpath = os.path.join(tmpdir, "graph")
        cmd = [
            "mpirun", "--oversubscribe", "-np", "1",
            binary, generator,
            "-s", str(seed),
            "-f", "edgelist",
            "-o", outpath,
            "-q",
        ] + args

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            raise RuntimeError(f"KaGen CLI failed: {result.stderr}")

        return parse_edgelist_file(outpath)


def parse_edgelist_file(path):
    """Parse KaGen DIMACS edgelist output file.

    Format:
        p <n> <m>
        e <u> <v>   (1-indexed)

    Returns 0-indexed sorted numpy array of shape (M, 2).
    """
    edges = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("e "):
                parts = line.split()
                u, v = int(parts[1]) - 1, int(parts[2]) - 1  # 1-indexed -> 0-indexed
                edges.append((u, v))

    if not edges:
        return np.empty((0, 2), dtype=np.uint64)

    arr = np.array(edges, dtype=np.uint64)
    idx = np.lexsort((arr[:, 1], arr[:, 0]))
    return arr[idx]


def sorted_edges(graph):
    """Return sorted edge array from a kagen.Graph."""
    edges = graph.edges()
    if edges.shape[0] == 0:
        return edges
    idx = np.lexsort((edges[:, 1], edges[:, 0]))
    return edges[idx]
