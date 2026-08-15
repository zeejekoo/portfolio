"""
Agent tool functions.

Each function here is a small, deterministic unit that answers ONE question
about the aging brain data. They are intentionally decoupled from the
LLM: they are just Python functions with typed args and structured
returns. `agent.py` wraps them with JSON Schema so Claude can call them.

Categories:
    Meta        - what's in the data?         (get_dataset_summary)
    Aggregate   - who is where, in numbers?   (get_cell_composition)
    Visualize   - show me a picture           (plot_spatial_map)
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Sequence

import matplotlib

matplotlib.use("Agg")  # headless server: never try to open a GUI window
import matplotlib.pyplot as plt
import numpy as np

from src.data_loader import Dataset, Age, load, subset, summary, PROJECT_ROOT


FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Meta - what's in the dataset?
# ---------------------------------------------------------------------------
def get_dataset_summary(dataset: Dataset) -> dict:
    """
    Return a compact summary of a dataset (cell count, ages, donors,
    cell type distribution). Used by the agent to orient itself before
    calling more expensive tools.
    """
    return summary(dataset)


# ---------------------------------------------------------------------------
# 2. Aggregate - cell type composition, filtered by conditions
# ---------------------------------------------------------------------------
def get_cell_composition(
    dataset: Dataset,
    age: Optional[Age] = None,
    donor_ids: Optional[Sequence[str]] = None,
    slices: Optional[Sequence[int]] = None,
    top_n: int = 15,
) -> dict:
    """
    Return cell type counts and percentages within the given subset.

    Any of age / donor_ids / slices can be omitted to skip that filter.
    The subset is described in the return value so agent responses can
    quote it back accurately.
    """
    adata = subset(dataset, age=age, donor_ids=donor_ids, slices=slices)

    counts = adata.obs["cell_type"].value_counts()
    total = int(counts.sum())
    pcts = (counts / total * 100).round(2)

    top = counts.head(top_n)
    top_pct = pcts.head(top_n)

    return {
        "dataset": dataset,
        "filters": {
            "age": age,
            "donor_ids": list(donor_ids) if donor_ids else None,
            "slices": list(slices) if slices else None,
        },
        "n_cells_in_subset": total,
        "n_cell_types": int((counts > 0).sum()),
        "composition_top_n": [
            {"cell_type": ct, "count": int(c), "percent": float(p)}
            for ct, c, p in zip(top.index, top.values, top_pct.values)
        ],
    }


# ---------------------------------------------------------------------------
# 3. Visualize - spatial map, colored by cell type
# ---------------------------------------------------------------------------
def plot_spatial_map(
    dataset: Dataset,
    slice_id: int,
    age: Optional[Age] = None,
    donor_id: Optional[str] = None,
    cell_types: Optional[Sequence[str]] = None,
    point_size: float = 0.5,
    output_name: Optional[str] = None,
) -> str:
    """
    Render a 2D spatial map of cells (colored by cell_type) for a given
    anatomical slice. Filters by age / donor / cell types if provided.

    Returns the path to the saved PNG (under figures/).

    Design note: we force one slice per call. Overlaying multiple slices
    in the same coordinate system produces an unreadable mess because
    slices are different physical brain sections.
    """
    adata = subset(
        dataset,
        age=age,
        donor_ids=[donor_id] if donor_id else None,
        cell_types=cell_types,
        slices=[slice_id],
    )

    if adata.n_obs == 0:
        raise ValueError(
            f"No cells match: dataset={dataset}, slice={slice_id}, "
            f"age={age}, donor={donor_id}, cell_types={cell_types}"
        )

    x = adata.obs["center_x"].values
    y = adata.obs["center_y"].values
    ct = adata.obs["cell_type"].astype(str).values

    # Deterministic color mapping so repeat calls look identical
    unique_ct = sorted(np.unique(ct))
    palette = plt.get_cmap("tab20", len(unique_ct))
    ct_to_color = {c: palette(i) for i, c in enumerate(unique_ct)}
    colors = [ct_to_color[c] for c in ct]

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(x, y, c=colors, s=point_size, linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()  # image conventions: y grows downward
    ax.set_xticks([])
    ax.set_yticks([])

    filters = [f"slice={slice_id}"]
    if age: filters.append(f"age={age}")
    if donor_id: filters.append(f"donor={donor_id}")
    if cell_types: filters.append(f"cell_types={len(cell_types)}")
    ax.set_title(f"{dataset} | " + ", ".join(filters), fontsize=10)

    # Compact legend on the right
    handles = [
        plt.Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=ct_to_color[c], markersize=6, label=c)
        for c in unique_ct
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
              fontsize=7, frameon=False)

    if output_name is None:
        stub = f"{dataset}_slice{slice_id}"
        if age: stub += f"_{age}"
        if donor_id: stub += f"_{donor_id.replace('MsBrainAgingSpatialDonor_','d')}"
        output_name = f"{stub}.png"

    out_path = FIGURES_DIR / output_name
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return str(out_path.relative_to(PROJECT_ROOT))


# ---------------------------------------------------------------------------
# Self-test - `python -m src.tools`
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import json

    print("=== 1) summary ===")
    print(json.dumps(get_dataset_summary("merfish_control"),
                     indent=2, default=str)[:400], "...\n")

    print("=== 2) composition (old, slice 0) ===")
    print(json.dumps(get_cell_composition("merfish_control", age="90wk",
                                          slices=[0], top_n=5),
                     indent=2, default=str))

    print("\n=== 3) spatial map (old, slice 0) ===")
    p = plot_spatial_map("merfish_control", slice_id=0, age="90wk")
    print(f"Saved: {p}")
