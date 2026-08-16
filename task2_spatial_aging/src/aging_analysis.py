"""
시험 Task 2 (2)-3: 나이별 특정 세포 타입 비교.

Young (4wk) vs Old (90wk) 마우스에서
- 세포 타입 비율 변화 (per-donor)
- 통계 유의성 (Mann-Whitney U)
- 공간 분포 시각화

Task 3 agent 가 tool 로 노출 예정.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from src.data_loader import Dataset, load, subset, PROJECT_ROOT


FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


def compare_celltype_ratio_by_age(
    dataset: Dataset = "merfish_control",
    cell_type: str = "microglial cell",
    young_age: str = "4wk",
    old_age: str = "90wk",
) -> dict:
    """
    두 나이 그룹 간 특정 세포 타입 비율 차이 정량 비교.

    각 donor 마다 (해당 세포 타입 셀 수 / 전체 셀 수) × 100 계산.
    donor 는 독립 표본이므로 Mann-Whitney U 로 그룹 간 차이 검정.

    반환:
        per_donor          : {age: {donor_id: percent}}
        young/old mean·std : 그룹 요약
        fold_change        : old/young 비
        p_value            : Mann-Whitney U 검정 (양측)
    """
    result = {
        "dataset": dataset, "cell_type": cell_type,
        "young_age": young_age, "old_age": old_age,
        "per_donor": {"young": {}, "old": {}},
    }

    for age_label, age in [("young", young_age), ("old", old_age)]:
        adata_age = subset(dataset, age=age)
        for donor in sorted(adata_age.obs["donor_id"].astype(str).unique()):
            donor_mask = adata_age.obs["donor_id"].astype(str) == donor
            donor_data = adata_age.obs[donor_mask]
            total = len(donor_data)
            ct_count = (donor_data["cell_type"].astype(str) == cell_type).sum()
            pct = ct_count / total * 100 if total > 0 else 0
            result["per_donor"][age_label][donor] = round(float(pct), 3)

    young_vals = list(result["per_donor"]["young"].values())
    old_vals = list(result["per_donor"]["old"].values())

    result["young_mean"] = float(np.mean(young_vals))
    result["young_std"] = float(np.std(young_vals, ddof=1))
    result["old_mean"] = float(np.mean(old_vals))
    result["old_std"] = float(np.std(old_vals, ddof=1))
    result["fold_change"] = round(result["old_mean"] / result["young_mean"], 3) \
        if result["young_mean"] > 0 else float("inf")

    stat, p = mannwhitneyu(young_vals, old_vals, alternative="two-sided")
    result["mw_stat"] = float(stat)
    result["p_value"] = float(p)

    return result


def plot_celltype_spatial_by_age(
    dataset: Dataset = "merfish_control",
    cell_type: str = "microglial cell",
    slice_id: int = 0,
    ages: Sequence[str] = ("4wk", "24wk", "90wk"),
    output_name: Optional[str] = None,
) -> str:
    """
    나이별 spatial map (1행 grid). 각 서브플롯:
    - 배경: 전체 세포 회색 흐리게
    - 강조: 지정 cell type 만 빨강

    각 나이에서 셀 수 가장 많은 donor 하나씩 대표로.
    """
    adata = load(dataset)
    slice_mask = adata.obs["slice"].astype(str) == str(slice_id)
    slice_adata = adata[slice_mask]

    fig, axes = plt.subplots(1, len(ages), figsize=(6 * len(ages), 6), squeeze=False)
    axes_flat = axes.flatten()

    for i, age in enumerate(ages):
        age_data = slice_adata[slice_adata.obs["age"].astype(str) == age]
        if age_data.n_obs == 0:
            axes_flat[i].set_title(f"{age} — 데이터 없음")
            axes_flat[i].axis("off")
            continue

        donor_counts = age_data.obs["donor_id"].value_counts()
        top_donor = str(donor_counts.index[0])
        donor_mask = age_data.obs["donor_id"].astype(str) == top_donor
        donor_data = age_data[donor_mask]

        x = donor_data.obs["center_x"].astype(float).values
        y = donor_data.obs["center_y"].astype(float).values
        ct = donor_data.obs["cell_type"].astype(str).values

        is_target = ct == cell_type

        ax = axes_flat[i]
        ax.scatter(x[~is_target], y[~is_target],
                   c="lightgray", s=1, alpha=0.3, linewidths=0)
        ax.scatter(x[is_target], y[is_target],
                   c="#d62728", s=8, alpha=0.9, linewidths=0)

        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_xticks([]); ax.set_yticks([])

        n_target = int(is_target.sum())
        n_total = len(ct)
        pct = n_target / n_total * 100 if n_total > 0 else 0
        short = top_donor.replace("MsBrainAgingSpatialDonor_", "d")
        ax.set_title(f"{age} — {short}\n"
                     f"{cell_type}: {n_target:,}/{n_total:,} ({pct:.2f}%)",
                     fontsize=11)

    fig.suptitle(f"{dataset} slice {slice_id} — "
                 f"{cell_type} spatial distribution by age",
                 fontsize=13, y=1.02)
    plt.tight_layout()

    if output_name is None:
        clean_ct = cell_type.replace(" ", "_")
        output_name = f"aging_{clean_ct}_slice{slice_id}.png"
    out_path = FIGURES_DIR / output_name
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return str(out_path.relative_to(PROJECT_ROOT))


if __name__ == "__main__":
    print("=" * 60)
    print("### 시험 (2)-3: Microglia 나이별 비교 (control) ###")
    print("=" * 60)
    result = compare_celltype_ratio_by_age(
        dataset="merfish_control",
        cell_type="microglial cell",
        young_age="4wk",
        old_age="90wk",
    )

    print(f"\n세포 타입: {result['cell_type']}")
    print(f"\n[Young ({result['young_age']})]  donor별 microglia 비율")
    for d, p in result["per_donor"]["young"].items():
        print(f"  {d}: {p:.3f}%")
    print(f"  → 평균 {result['young_mean']:.3f}% ± {result['young_std']:.3f}%")

    print(f"\n[Old ({result['old_age']})]  donor별 microglia 비율")
    for d, p in result["per_donor"]["old"].items():
        print(f"  {d}: {p:.3f}%")
    print(f"  → 평균 {result['old_mean']:.3f}% ± {result['old_std']:.3f}%")

    print(f"\n=== 통계 검정 ===")
    print(f"  Fold change (old/young): {result['fold_change']}x")
    print(f"  Mann-Whitney U p-value : {result['p_value']:.4f}")
    if result["p_value"] < 0.05:
        print(f"  → 유의미한 차이 (p<0.05)")
    else:
        print(f"  → 통계적 유의성 미확보 (p>=0.05)  ※ n=4 vs 4 소표본")

    print(f"\n=== 공간 지도 생성 ===")
    for sl in [0, 1]:
        p = plot_celltype_spatial_by_age(
            dataset="merfish_control",
            cell_type="microglial cell",
            slice_id=sl,
        )
        print(f"  slice {sl}: {p}")
