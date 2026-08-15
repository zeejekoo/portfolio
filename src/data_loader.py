"""
Data loading utilities for the aging brain MERFISH + snRNA-seq demo.

All h5ad access flows through this module so paths, caching, and
subset conventions live in one place. Agent tools, notebooks, and
the Streamlit app should import from here instead of calling
scanpy.read_h5ad directly.

Data schema (confirmed via sanity check):
    obs: age {'4wk','24wk','90wk'}, donor_id (12 mice), slice {0,1,2},
         cell_type (13 broad types), center_x/center_y
    obsm: 'spatial_coords', 'X_umap', 'X_pca'
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Literal, Optional, Sequence

import anndata as ad
import scanpy as sc


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILES = {
    "merfish_control": PROJECT_ROOT / "merfish_control.h5ad",
    "merfish_lps":     PROJECT_ROOT / "merfish_lps.h5ad",
    "snrnaseq":        PROJECT_ROOT / "snrnaseq.h5ad",
}

Dataset = Literal["merfish_control", "merfish_lps", "snrnaseq"]
Age = Literal["4wk", "24wk", "90wk"]


@lru_cache(maxsize=3)
def load(dataset: Dataset) -> ad.AnnData:
    """Load one of the h5ad files, cached in memory for the process."""
    path = DATA_FILES[dataset]
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run scripts/download_data.sh first."
        )
    return sc.read_h5ad(path)


def subset(
    dataset: Dataset,
    age: Optional[Age | Sequence[Age]] = None,
    cell_types: Optional[Sequence[str]] = None,
    donor_ids: Optional[Sequence[str]] = None,
    slices: Optional[Sequence[int]] = None,
) -> ad.AnnData:
    """
    Return a fresh AnnData restricted to the given age(s), cell type(s),
    donor(s), and/or anatomical slice(s). None = no filter on that axis.

    The returned object is always a copy - safe to mutate.
    """
    adata = load(dataset)
    mask = None

    def _combine(m, new):
        return new if m is None else (m & new)

    if age is not None:
        ages = [age] if isinstance(age, str) else list(age)
        mask = _combine(mask, adata.obs["age"].isin(ages).values)
    if cell_types is not None:
        mask = _combine(mask, adata.obs["cell_type"].isin(cell_types).values)
    if donor_ids is not None:
        mask = _combine(mask, adata.obs["donor_id"].isin(donor_ids).values)
    if slices is not None:
        mask = _combine(mask, adata.obs["slice"].isin(slices).values)

    if mask is None:
        return adata.copy()
    return adata[mask].copy()


def summary(dataset: Dataset) -> dict:
    """Compact summary dict for logging / status bar / agent context."""
    adata = load(dataset)
    return {
        "dataset": dataset,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "ages": sorted(adata.obs["age"].unique().tolist()),
        "n_donors": int(adata.obs["donor_id"].nunique()),
        "n_slices": int(adata.obs["slice"].nunique()) if "slice" in adata.obs else None,
        "cell_types": adata.obs["cell_type"].value_counts().to_dict(),
    }


if __name__ == "__main__":
    import json
    s = summary("merfish_control")
    s["cell_types"] = dict(list(s["cell_types"].items())[:5]) | {"...": "..."}
    print(json.dumps(s, indent=2, default=str))
