#!/usr/bin/env bash
# Thin wrapper around `python -m src.pipeline fetch all`.
#
# Activates the local venv if it exists, otherwise falls back to system python3.
# Forwards every argument to the pipeline CLI, so you can still pass --overwrite,
# --limit, etc.
#
# Usage:
#   ./scripts/fetch_all.sh                # fetch everything
#   ./scripts/fetch_all.sh --limit 3      # smoke test: 3 MAGs per BioProject
#   ./scripts/fetch_all.sh --overwrite    # ignore existing artifacts

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}"

if [[ -x ".venv/bin/python" ]]; then
  PYTHON=".venv/bin/python"
else
  PYTHON="python3"
fi

if [[ -z "${NCBI_EMAIL:-}" ]]; then
  echo "[fetch_all.sh] WARNING: NCBI_EMAIL is unset. Set it before large NCBI batches." >&2
fi

exec "${PYTHON}" -m src.pipeline fetch all "$@"
