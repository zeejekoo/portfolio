"""
Agent tool 함수 모음.

카테고리:
    Meta      - get_dataset_summary
    Aggregate - get_cell_composition
    Visualize - plot_spatial_map, plot_all_samples_grid (color_by 지원)
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from src.data_loader import Dataset, Age, load, subset, summary, PROJECT_ROOT


FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# 팔레트
# ---------------------------------------------------------------------------
VIVID_PALETTE = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
    "#9a6324", "#800000", "#000075", "#808000", "#ffe119",
]


def _sort_categories(cats) -> list:
    """모두 정수 파싱되면 숫자 순, 아니면 사전 순."""
    cats = [str(c) for c in cats]
    try:
        return sorted(cats, key=lambda x: int(x))
    except ValueError:
        return sorted(cats)


def _color_map(categories) -> dict:
    """카테고리 → 색 결정적 매핑. 15 이하는 VIVID_PALETTE, 초과 시 tab20/tab20b 확장."""
    cats = _sort_categories(categories)
    n = len(cats)
    if n <= len(VIVID_PALETTE):
        return {c: VIVID_PALETTE[i] for i, c in enumerate(cats)}
    # 확장 팔레트: 15 + tab20(20) + tab20b(20) = 55색
    extra = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors)
    all_colors = list(VIVID_PALETTE) + extra
    return {c: all_colors[i % len(all_colors)] for i, c in enumerate(cats)}


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
    point_size: float = 2.0,
    output_name: Optional[str] = None,
) -> str:
    adata = subset(
        dataset, age=age,
        donor_ids=[donor_id] if donor_id else None,
        cell_types=cell_types, slices=[slice_id],
    )
    if adata.n_obs == 0:
        raise ValueError(f"조건에 맞는 세포 없음: {dataset}, slice={slice_id}")

    x = adata.obs["center_x"].values
    y = adata.obs["center_y"].values
    ct = adata.obs["cell_type"].astype(str).values

    all_ct = load(dataset).obs["cell_type"].astype(str).unique()
    ct_to_color = _color_map(all_ct)
    colors = [ct_to_color[c] for c in ct]

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(x, y, c=colors, s=point_size, linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_xticks([]); ax.set_yticks([])

    filters = [f"slice={slice_id}"]
    if age: filters.append(f"age={age}")
    if donor_id: filters.append(f"donor={donor_id}")
    ax.set_title(f"{dataset} | " + ", ".join(filters), fontsize=10)

    unique_ct = _sort_categories(set(ct))
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=ct_to_color[c], markersize=8, label=c)
               for c in unique_ct]
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
# 4. Visualize — 한 슬라이스 모든 donor grid (color_by 지원)
# ---------------------------------------------------------------------------
def plot_all_samples_grid(
    dataset: Dataset,
    slice_id: int,
    color_by: str = "cell_type",   # "cell_type" or "leiden"
    resolution: float = 1.0,       # leiden 캐시 파일 특정용
    cols: int = 4,
    point_size: float = 2.5,
    output_name: Optional[str] = None,
) -> str:
    """
    한 슬라이스 안의 모든 donor 를 grid 로 그림.
    color_by="cell_type" 이면 저자 라벨 색상, "leiden" 이면 캐시된 클러스터 색상.
    """
    # 데이터 소스: leiden 이면 캐시에서 로드, cell_type 이면 원본
    if color_by == "leiden":
        from src.clustering import _cache_path
        cache = _cache_path(dataset, resolution=resolution)
        if not cache.exists():
            raise FileNotFoundError(
                f"Leiden 캐시 없음: {cache.name}\n"
                f"먼저 실행: python -m src.clustering"
            )
        adata_all = sc.read_h5ad(cache)
    else:
        adata_all = load(dataset)

    if color_by not in adata_all.obs.columns:
        raise KeyError(f"obs 에 '{color_by}' 없음")

    for col in ["slice", "donor_id", "center_x", "center_y"]:
        if col not in adata_all.obs.columns:
            raise KeyError(f"obs 에 '{col}' 없음")

    has_age = "age" in adata_all.obs.columns

    slice_mask = adata_all.obs["slice"].astype(str) == str(slice_id)
    if slice_mask.sum() == 0:
        raise ValueError(f"{dataset} 에 slice {slice_id} 없음")
    slice_data = adata_all[slice_mask]

    donors_meta = (
        slice_data.obs[["donor_id"] + (["age"] if has_age else [])]
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

    all_cats = _sort_categories(adata_all.obs[color_by].astype(str).unique())
    cat_to_color = _color_map(all_cats)

    fig, axes = plt.subplots(n_rows, cols, figsize=(cols * 4, n_rows * 4), squeeze=False)
    axes_flat = axes.flatten()

    for i, row in enumerate(donors_meta.itertuples(index=False)):
        donor_id = row.donor_id
        age = getattr(row, "age", None)
        ax = axes_flat[i]

        donor_mask = slice_data.obs["donor_id"].astype(str) == str(donor_id)
        sub = slice_data[donor_mask]

        x = sub.obs["center_x"].values
        y = sub.obs["center_y"].values
        cat = sub.obs[color_by].astype(str).values
        colors = [cat_to_color[c] for c in cat]

        ax.scatter(x, y, c=colors, s=point_size, linewidths=0)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])

        short = donor_id.replace("MsBrainAgingSpatialDonor_", "d")
        title = f"{short} ({age})" if age else short
        ax.set_title(title, fontsize=10)

    for j in range(n_donors, len(axes_flat)):
        axes_flat[j].axis("off")

    # legend (많으면 2열)
    ncol_legend = 2 if len(all_cats) > 15 else 1
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=cat_to_color[c], markersize=8, label=c)
               for c in all_cats]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=8 if ncol_legend == 2 else 9, frameon=False, ncol=ncol_legend)

    fig.suptitle(f"{dataset} — slice {slice_id} — colored by {color_by} — {n_donors} donors",
                 fontsize=14, y=1.01)
    plt.tight_layout()

    if output_name is None:
        output_name = f"{dataset}_all_donors_slice{slice_id}_{color_by}.png"
    out_path = FIGURES_DIR / output_name
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(out_path.relative_to(PROJECT_ROOT))


if __name__ == "__main__":
    print("=== leiden 지도 self-test (control slice 0) ===")
    p = plot_all_samples_grid("merfish_control", slice_id=0, color_by="leiden")
    print(f"저장: {p}")
