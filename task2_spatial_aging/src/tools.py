"""
Agent tool 함수 모음.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.data_loader import Dataset, Age, load, subset, summary, PROJECT_ROOT


FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# vivid 팔레트 (tab20 은 파스텔 섞여있어 옅게 보임)
# ---------------------------------------------------------------------------
VIVID_PALETTE = [
    "#e6194b",  # red
    "#3cb44b",  # green
    "#4363d8",  # blue
    "#f58231",  # orange
    "#911eb4",  # purple
    "#42d4f4",  # cyan
    "#f032e6",  # magenta
    "#bfef45",  # lime
    "#fabed4",  # pink
    "#469990",  # teal
    "#9a6324",  # brown
    "#800000",  # maroon
    "#000075",  # navy
    "#808000",  # olive
    "#ffe119",  # yellow (마지막에 배치 — 흰 배경 대비 애매)
]


def _color_map(categories: Sequence[str]) -> dict:
    """세포 타입 → 색 결정적 매핑 (반복 호출해도 같은 색)."""
    cats = sorted(set(categories))
    return {c: VIVID_PALETTE[i % len(VIVID_PALETTE)] for i, c in enumerate(cats)}


# ---------------------------------------------------------------------------
# 1. Meta
# ---------------------------------------------------------------------------
def get_dataset_summary(dataset: Dataset) -> dict:
    return summary(dataset)


# ---------------------------------------------------------------------------
# 2. Aggregate
# ---------------------------------------------------------------------------
def get_cell_composition(
    dataset: Dataset,
    age: Optional[Age] = None,
    donor_ids: Optional[Sequence[str]] = None,
    slices: Optional[Sequence[int]] = None,
    top_n: int = 15,
) -> dict:
    adata = subset(dataset, age=age, donor_ids=donor_ids, slices=slices)
    counts = adata.obs["cell_type"].value_counts()
    total = int(counts.sum())
    pcts = (counts / total * 100).round(2) if total > 0 else counts * 0
    top = counts.head(top_n)
    top_pct = pcts.head(top_n)
    return {
        "dataset": dataset,
        "filters": {"age": age, "donor_ids": list(donor_ids) if donor_ids else None,
                    "slices": list(slices) if slices else None},
        "n_cells_in_subset": total,
        "n_cell_types": int((counts > 0).sum()),
        "composition_top_n": [
            {"cell_type": ct, "count": int(c), "percent": float(p)}
            for ct, c, p in zip(top.index, top.values, top_pct.values)
        ],
    }


# ---------------------------------------------------------------------------
# 3. Visualize — 단일 (donor, slice)
# ---------------------------------------------------------------------------
def plot_spatial_map(
    dataset: Dataset,
    slice_id: int,
    age: Optional[Age] = None,
    donor_id: Optional[str] = None,
    cell_types: Optional[Sequence[str]] = None,
    point_size: float = 2.0,  # 이전 0.5 → 2.0
    output_name: Optional[str] = None,
) -> str:
    adata = subset(
        dataset, age=age,
        donor_ids=[donor_id] if donor_id else None,
        cell_types=cell_types, slices=[slice_id],
    )
    if adata.n_obs == 0:
        raise ValueError(
            f"조건에 맞는 세포가 없음: dataset={dataset}, slice={slice_id}"
        )

    x = adata.obs["center_x"].values
    y = adata.obs["center_y"].values
    ct = adata.obs["cell_type"].astype(str).values

    # dataset 전체 cell_type 기준 색 매핑 (subplot 간 일관)
    all_ct = load(dataset).obs["cell_type"].astype(str).unique()
    ct_to_color = _color_map(all_ct)
    colors = [ct_to_color[c] for c in ct]

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(x, y, c=colors, s=point_size, linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])

    filters = [f"slice={slice_id}"]
    if age: filters.append(f"age={age}")
    if donor_id: filters.append(f"donor={donor_id}")
    ax.set_title(f"{dataset} | " + ", ".join(filters), fontsize=10)

    unique_ct_in_subset = sorted(set(ct))
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=ct_to_color[c], markersize=8, label=c)
               for c in unique_ct_in_subset]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
              fontsize=8, frameon=False)

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
# 4. Visualize — 한 슬라이스, 모든 donor grid (시험 (1)-1)
# ---------------------------------------------------------------------------
def plot_all_samples_grid(
    dataset: Dataset,
    slice_id: int,
    cols: int = 4,
    point_size: float = 1.5,  # 이전 0.3 → 1.5
    output_name: Optional[str] = None,
) -> str:
    """
    한 슬라이스 안의 모든 donor 를 grid 로 그린 대형 PNG.
    나이 순 정렬(4wk → 24wk → 90wk), vivid 팔레트, subplot 간 색 일관.
    """
    adata_all = load(dataset)

    for col in ["slice", "donor_id", "cell_type", "center_x", "center_y"]:
        if col not in adata_all.obs.columns:
            raise KeyError(f"{dataset} 의 obs 에 '{col}' 컬럼 없음")

    has_age = "age" in adata_all.obs.columns
    slice_mask = adata_all.obs["slice"].astype(str) == str(slice_id)
    if slice_mask.sum() == 0:
        raise ValueError(f"{dataset} 에 slice {slice_id} 데이터 없음")

    donors_meta = (
        adata_all.obs.loc[slice_mask, ["donor_id"] + (["age"] if has_age else [])]
        .drop_duplicates()
    )
    if has_age:
        donors_meta = donors_meta.assign(
            _age_num=donors_meta["age"].astype(str).str.replace("wk", "").astype(int)
        ).sort_values(["_age_num", "donor_id"]).drop(columns=["_age_num"])
    else:
        donors_meta = donors_meta.sort_values("donor_id")

    n_donors = len(donors_meta)
    n_rows = (n_donors + cols - 1) // cols

    # 통일 팔레트 (dataset 전체 cell_type 기준)
    ct_to_color = _color_map(adata_all.obs["cell_type"].astype(str).unique())

    fig, axes = plt.subplots(n_rows, cols, figsize=(cols * 4, n_rows * 4), squeeze=False)
    axes_flat = axes.flatten()

    for i, row in enumerate(donors_meta.itertuples(index=False)):
        donor_id = row.donor_id
        age = getattr(row, "age", None)
        ax = axes_flat[i]

        sub = subset(dataset, donor_ids=[donor_id], slices=[slice_id])
        x = sub.obs["center_x"].values
        y = sub.obs["center_y"].values
        ct = sub.obs["cell_type"].astype(str).values
        colors = [ct_to_color[c] for c in ct]

        ax.scatter(x, y, c=colors, s=point_size, linewidths=0)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])

        short = donor_id.replace("MsBrainAgingSpatialDonor_", "d")
        title = f"{short} ({age})" if age else short
        ax.set_title(title, fontsize=10)

    for j in range(n_donors, len(axes_flat)):
        axes_flat[j].axis("off")

    # 전체 legend
    all_ct = sorted(adata_all.obs["cell_type"].astype(str).unique())
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=ct_to_color[c], markersize=9, label=c)
               for c in all_ct]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=9, frameon=False)

    fig.suptitle(f"{dataset} — slice {slice_id} — {n_donors} donors",
                 fontsize=14, y=1.01)
    plt.tight_layout()

    if output_name is None:
        output_name = f"{dataset}_all_donors_slice{slice_id}.png"
    out_path = FIGURES_DIR / output_name
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)

    return str(out_path.relative_to(PROJECT_ROOT))


# ---------------------------------------------------------------------------
# self-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    p = plot_all_samples_grid("merfish_lps", slice_id=1)
    print(f"저장: {p}")
