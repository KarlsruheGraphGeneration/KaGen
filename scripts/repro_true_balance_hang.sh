#!/usr/bin/env bash
#SBATCH --nodes=16
#SBATCH -o repro_true_balance_hang-%j-log.txt
#SBATCH -e repro_true_balance_hang-%j-err.txt
#SBATCH -J kagen-true-balance-hang-repro
#SBATCH --partition=micro
#SBATCH --time=0-00:15:00
#SBATCH --get-user-env
#SBATCH --account=pn72pu
#SBATCH --switches=1
#SBATCH --constraint='work'
#SBATCH --ear=off
#
# Isolates the KaCCv2 hang to KaGen's graph generation/redistribution step, without
# involving KaCCv2 or the MultiStep algorithm at all. Submit with `sbatch` on
# SuperMUC (same partition/account/module setup as the original KaCCv2 job); the
# #SBATCH directives above and the module/env setup below mirror
#   KaCCv2/eval/data/supermuc/multistep-reachability-comm-modes-rmat-balanced_26_07_22/jobfiles/
#     multistep-reachability-comm-modes-rmat-balanced-in2_rmat_n18-p768
# which is the job that hung.
#
# If this hangs too: the bug is in KaGen (most likely in
# RedistributeEdgesTrueBalance(), kagen/tools/postprocessor.cpp), not in KaCCv2.
# That function now prints root-rank "breadcrumb" progress lines to stderr; the last
# one printed without a matching "after ..." pinpoints which MPI collective never
# returned -- see repro_true_balance_hang-<jobid>-log.txt / -err.txt (or the console,
# if run interactively).
#
# Usage:
#   sbatch scripts/repro_true_balance_hang.sh
#   KAGEN_BINARY=/path/to/build/app/KaGen sbatch scripts/repro_true_balance_hang.sh
#   NTASKS=48 PERHOST=48 N=1000000 M=8000000 sbatch scripts/repro_true_balance_hang.sh   # small smoke test
#
# To run interactively inside an existing salloc instead of via sbatch, just
# `bash scripts/repro_true_balance_hang.sh` after exporting the same variables --
# the #SBATCH lines are ignored outside of sbatch.
#
# Requires a KaGen CLI binary already built with -DKAGEN_BUILD_APPS=ON from a
# checkout that includes the breadcrumb instrumentation (this repo). By default
# this resolves to build/app/KaGen next to this checkout's own scripts/ directory
# (i.e. wherever you cloned KaGen on this system); override KAGEN_BINARY if your
# build directory is named/located differently. Build it with, e.g.:
#   cmake -S . -B build -DKAGEN_BUILD_APPS=ON -DKAGEN_USE_XXHASH=ON \
#         -DKAGEN_USE_CGAL=OFF -DKAGEN_USE_SPARSEHASH=OFF
#   cmake --build build --target KaGenApp -j8

set -euo pipefail

# --- Configuration (override via environment) --------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KAGEN_BINARY="${KAGEN_BINARY:-$SCRIPT_DIR/../build/app/KaGen}"

# Same option string as the hanging KaCCv2 run.
OPTION_STRING="${OPTION_STRING:-type=rmat;a=0.57;b=0.19;c=0.19;redistribution=balance-edges-strict;n=201326592;m=1610612736}"

NTASKS="${NTASKS:-768}"
PERHOST="${PERHOST:-48}"

if [[ ! -x "$KAGEN_BINARY" ]]; then
    echo "[repro] error: KAGEN_BINARY '$KAGEN_BINARY' not found or not executable." >&2
    echo "[repro] Build it first (see the comment block at the top of this script) or set KAGEN_BINARY." >&2
    exit 1
fi

module restore kaccv2
module load slurm_setup

echo "[repro] Binary:  $KAGEN_BINARY"
echo "[repro] Ranks:   $NTASKS (--perhost $PERHOST)"
echo "[repro] Options: $OPTION_STRING"
echo "[repro] Starting at $(date)"

export I_MPI_PIN_CELL=core
export I_MPI_PIN_DOMAIN=1:compact
export I_MPI_JOB_TIMEOUT="${I_MPI_JOB_TIMEOUT:-120}"

mpiexec -n "$NTASKS" --perhost "$PERHOST" \
    "$KAGEN_BINARY" options "$OPTION_STRING" \
    --stats basic \
    -f noop

echo "[repro] Finished at $(date)"
