"""
Agent tool 함수 모음.

카테고리:
    Meta        - get_dataset_summary
    Aggregate   - get_cell_composition
    Visualize   - plot_spatial_map (단일 샘플)
                  plot_all_samples_grid (한 슬라이스, 모든 donor)
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
# 1. Meta
# ---------------------------------------------------------------------------
def get_dataset_summary(dataset: Dataset) -> dict:
    """agent 가 데이터셋 방향을 잡을 때 첫 호출용 압축 요약."""
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
    """조건별 세포 타입 조성. 필터 조합 자유."""
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
    point_size: float = 0.5,
    output_name: Optional[str] = None,
) -> str:
    """단일 슬라이스의 2D spatial map (cell_type 색상)."""
    adata = subset(
        dataset, age=age,
        donor_ids=[donor_id] if donor_id else None,
        cell_types=cell_types, slices=[slice_id],
    )
    if adata.n_obs == 0:
        raise ValueError(
            f"조건에 맞는 세포가 없음: dataset={dataset}, slice={slice_id}, "
            f"age={age}, donor={donor_id}"
        )

    x = adata.obs["center_x"].values
    y = adata.obs["center_y"].values
    ct = adata.obs["cell_type"].astype(str).values

    unique_ct = sorted(np.unique(ct))
    palette = plt.get_cmap("tab20", len(unique_ct))
    ct_to_color = {c: palette(i) for i, c in enumerate(unique_ct)}
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

    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=ct_to_color[c], markersize=6, label=c)
               for c in unique_ct]
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
# 4. Visualize — 한 슬라이스, 모든 donor grid (시험 (1)-1)
# ---------------------------------------------------------------------------
def plot_all_samples_grid(
    dataset: Dataset,
    slice_id: int,
    cols: int = 4,
    point_size: float = 0.3,
    output_name: Optional[str] = None,
) -> str:
    """
    한 슬라이스 안의 모든 donor 를 grid 로 그린 대형 PNG.
    나이 순으로 정렬(4wk → 24wk → 90wk)해 시각적 스토리 형성.

    파라미터:
        dataset      : merfish_control / merfish_lps
        slice_id     : 0, 1, 2 중 하나
        cols         : grid 열 수
        point_size   : 점 크기
        output_name  : 저장 파일명 (미지정시 자동)

    반환:
        저장된 PNG 의 프로젝트 루트 기준 상대 경로.
    """
    adata_all = load(dataset)

    # 필수 컬럼 존재 확인 (merfish_lps 스키마가 다를 수 있음)
    for col in ["slice", "donor_id", "cell_type", "center_x", "center_y"]:
        if col not in adata_all.obs.columns:
            raise KeyError(f"{dataset} 의 obs 에 '{col}' 컬럼 없음")

    # age 는 없을 수도 있음 (있으면 나이 순 정렬, 없으면 donor_id 순)
    has_age = "age" in adata_all.obs.columns

    slice_mask = adata_all.obs["slice"].astype(str) == str(slice_id)
    if slice_mask.sum() == 0:
        raise ValueError(f"{dataset} 에 slice {slice_id} 데이터 없음")

    donors_meta = (
        adata_all.obs.loc[slice_mask, ["donor_id"] + (["age"] if has_age else [])]
        .drop_duplicates()
    )
    if has_age:
        # 4wk → 24wk → 90wk 로 정렬 (숫자 파싱)
        donors_meta = donors_meta.assign(
            _age_num=donors_meta["age"].astype(str).str.replace("wk", "").astype(int)
        ).sort_values(["_age_num", "donor_id"]).drop(columns=["_age_num"])
    else:
        donors_meta = donors_meta.sort_values("donor_id")

    n_donors = len(donors_meta)
    n_rows = (n_donors + cols - 1) // cols

    # 전역 색 팔레트 (dataset 전체 cell_type 기준 → 서브플롯 간 일관됨)
    unique_ct = sorted(adata_all.obs["cell_type"].astype(str).unique())
    palette = plt.get_cmap("tab20", len(unique_ct))
    ct_to_color = {c: palette(i) for i, c in enumerate(unique_ct)}

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

    # 남은 빈 subplot 숨김
    for j in range(n_donors, len(axes_flat)):
        axes_flat[j].axis("off")

    # 전체 legend (fig 오른쪽 바깥)
    handles = [plt.Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=ct_to_color[c], markersize=8, label=c)
               for c in unique_ct]
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
    import json

    print("=== 1) summary ===")
    print(json.dumps(get_dataset_summary("merfish_control"),
                     indent=2, default=str, ensure_ascii=False)[:300], "...\n")

    print("=== 2) composition (old, slice 0) ===")
    print(json.dumps(get_cell_composition("merfish_control", age="90wk",
                                          slices=[0], top_n=5),
                     indent=2, default=str, ensure_ascii=False))

    print("\n=== 3) plot_spatial_map (single sample) ===")
    p1 = plot_spatial_map("merfish_control", slice_id=0, age="90wk")
    print(f"저장: {p1}")

    print("\n=== 4) plot_all_samples_grid (slice 0, all donors) ===")
    p2 = plot_all_samples_grid("merfish_control", slice_id=0)
    print(f"저장: {p2}")
