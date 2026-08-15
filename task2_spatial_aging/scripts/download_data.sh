#!/usr/bin/env bash
# Download MERFISH + snRNA-seq data for the aging brain demo.
# Source: CELLxGENE collection 31937775-0602-4e52-a799-b6acdd2bac2e
# Paper: Allen et al., Cell 2023, doi:10.1016/j.cell.2022.12.010
#
# Usage: bash scripts/download_data.sh
# Re-run safe: wget -c resumes/skips already-complete files.

set -euo pipefail

DATA_DIR="$(dirname "$0")/.."
cd "$DATA_DIR"

echo "[1/3] MERFISH control (351 MB)"
wget -c -O merfish_control.h5ad \
  https://datasets.cellxgene.cziscience.com/c93d78c2-ee17-4504-8d1c-17cf093ad7b5.h5ad

echo "[2/3] MERFISH LPS (412 MB)"
wget -c -O merfish_lps.h5ad \
  https://datasets.cellxgene.cziscience.com/9a610e48-3fb3-4813-9cd0-ec5ae190e49c.h5ad

echo "[3/3] snRNA-seq (2.53 GB)"
wget -c -O snrnaseq.h5ad \
  https://datasets.cellxgene.cziscience.com/752ed29d-bbb0-45bd-9c57-19edf68ef73e.h5ad

echo "Done. Files:"
ls -lh *.h5ad
