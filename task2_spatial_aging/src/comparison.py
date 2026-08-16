"""
(3)-1: 노화 vs LPS 공통 DEG.

노화 DEG (control 24wk vs 90wk) 와
LPS  DEG (control 24wk vs lps 24wk) 의 교집합.

핵심 설계:
    - adata.raw (원본 카운트) 에서 재정규화 (저자 스케일링된 X 사용 X)
    - Wilcoxon rank-sum test (scanpy rank_genes_groups)
    - 유의성 기준: pvals_adj < 0.05 AND |log2FC| > 0.5
    - 교집합의 통계적 유의성 = hypergeometric test
"""

from __future__ import annotations
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import hypergeom

from src.data_loader import Dataset, load, subset, PROJECT_ROOT


OUTPUTS_DIR = PROJECT_ROOT / "outputs"
OUTPUTS_DIR.mkdir(exist_ok=True)


def _prepare_for_deg(adata: ad.AnnData) -> ad.AnnData:
    """
    adata.raw (원본 카운트) 에서 다시 정규화+log 하여 DEG 계산 대비.
    저자 X (스케일링됨) 사용 시 log fold change 왜곡 방지.
    """
    if adata.raw is None:
        raise ValueError("adata.raw 없음 — raw counts 필요")

    # raw 슬롯의 counts + var 를 새 AnnData 로
    new_adata = ad.AnnData(
        X=adata.raw.X.copy(),
        obs=adata.obs.copy(),
        var=adata.raw.var.copy(),
    )
    sc.pp.normalize_total(new_adata, target_sum=1e4)
    sc.pp.log1p(new_adata)
    return new_adata


def _extract_deg(
    adata: ad.AnnData,
    group_col: str,
    target_group: str,
    reference_group: str,
    pval_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
) -> dict:
    """
    Wilcoxon rank-sum DEG. target vs reference.
    반환: up (target 에서 상승), down (target 에서 감소), 전체 결과 DataFrame.
    """
    sc.tl.rank_genes_groups(
        adata, group_col, groups=[target_group], reference=reference_group,
        method="wilcoxon", pts=True,
    )
    key = "rank_genes_groups"
    var_names = adata.var.get("feature_name", adata.var_names).astype(str).values

    # var_names_map: index → feature_name
    name_map = {i: n for i, n in zip(adata.var_names, var_names)}

    result = pd.DataFrame({
        "gene_id": adata.uns[key]["names"][target_group],
        "logfoldchange": adata.uns[key]["logfoldchanges"][target_group],
        "pvals_adj": adata.uns[key]["pvals_adj"][target_group],
    })
    result["gene"] = result["gene_id"].map(name_map).fillna(result["gene_id"])

    sig = result[
        (result["pvals_adj"] < pval_threshold) &
        (result["logfoldchange"].abs() > lfc_threshold)
    ]
    up = sig[sig["logfoldchange"] > 0]["gene"].tolist()
    down = sig[sig["logfoldchange"] < 0]["gene"].tolist()

    return {
        "up": up, "down": down,
        "n_up": len(up), "n_down": len(down),
        "all_result": result,
    }


def find_aging_deg(
    dataset: Dataset = "merfish_control",
    young_age: str = "24wk",
    old_age: str = "90wk",
    pval_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
) -> dict:
    """노화 DEG: young 대비 old 에서 up/down."""
    adata_raw = subset(dataset, age=[young_age, old_age])
    adata = _prepare_for_deg(adata_raw)
    adata.obs["_age"] = adata.obs["age"].astype(str)
    return _extract_deg(adata, "_age", old_age, young_age,
                        pval_threshold, lfc_threshold)


def find_lps_deg(
    age: str = "24wk",
    pval_threshold: float = 0.05,
    lfc_threshold: float = 0.5,
) -> dict:
    """LPS DEG: control 대비 lps 에서 up/down (age-matched)."""
    c = subset("merfish_control", age=age)
    l = subset("merfish_lps", age=age)

    c_prep = _prepare_for_deg(c)
    l_prep = _prepare_for_deg(l)
    c_prep.obs["_cond"] = "control"
    l_prep.obs["_cond"] = "lps"

    # 두 데이터셋 공통 유전자만 (feature_name 기준)
    common_genes = sorted(
        set(c_prep.var["feature_name"].astype(str)) &
        set(l_prep.var["feature_name"].astype(str))
    )
    c_prep = c_prep[:, c_prep.var["feature_name"].astype(str).isin(common_genes)]
    l_prep = l_prep[:, l_prep.var["feature_name"].astype(str).isin(common_genes)]

    merged = ad.concat([c_prep, l_prep], axis=0, join="inner")
    merged.var = c_prep.var.loc[merged.var_names].copy()
    return _extract_deg(merged, "_cond", "lps", "control",
                        pval_threshold, lfc_threshold)


def find_common_deg(
    direction: str = "any",
    aging_kwargs: Optional[dict] = None,
    lps_kwargs: Optional[dict] = None,
) -> dict:
    """
    노화 DEG ∩ LPS DEG.
    direction: 'up' (양쪽 다 상승), 'down' (양쪽 다 감소), 'any' (방향 무관 겹침)
    """
    aging = find_aging_deg(**(aging_kwargs or {}))
    lps = find_lps_deg(**(lps_kwargs or {}))

    if direction == "up":
        aging_set = set(aging["up"])
        lps_set = set(lps["up"])
    elif direction == "down":
        aging_set = set(aging["down"])
        lps_set = set(lps["down"])
    elif direction == "any":
        aging_set = set(aging["up"]) | set(aging["down"])
        lps_set = set(lps["up"]) | set(lps["down"])
    else:
        raise ValueError(f"direction must be up/down/any, got {direction}")

    common = aging_set & lps_set

    # Hypergeom test — universe = MERFISH panel 전체 (374)
    M = 374
    n = len(aging_set)
    N = len(lps_set)
    k = len(common)
    expected = N * n / M if M > 0 else 0
    p_value = float(hypergeom.sf(k - 1, M, n, N)) if k > 0 else 1.0
    enrichment = k / expected if expected > 0 else float("inf")

    return {
        "direction": direction,
        "n_aging_deg": len(aging_set),
        "n_lps_deg": len(lps_set),
        "n_common": len(common),
        "common_genes": sorted(common),
        "aging_only": sorted(aging_set - lps_set),
        "lps_only": sorted(lps_set - aging_set),
        "expected_random": round(expected, 2),
        "enrichment": round(enrichment, 2),
        "p_value": p_value,
        "aging_up": aging["up"],
        "aging_down": aging["down"],
        "lps_up": lps["up"],
        "lps_down": lps["down"],
    }


if __name__ == "__main__":
    print("=" * 60)
    print("### (3)-1: 노화 vs LPS 공통 발현 변화 유전자 ###")
    print("=" * 60)

    print("\n[노화 DEG 계산: control 24wk vs 90wk] ...")
    aging = find_aging_deg()
    print(f"  ↑ up  : {aging['n_up']}개 유전자")
    print(f"  ↓ down: {aging['n_down']}개")
    print(f"  예시 up   : {aging['up'][:10]}")
    print(f"  예시 down : {aging['down'][:10]}")

    print("\n[LPS DEG 계산: control 24wk vs lps 24wk] ...")
    lps = find_lps_deg()
    print(f"  ↑ up  : {lps['n_up']}개 유전자")
    print(f"  ↓ down: {lps['n_down']}개")
    print(f"  예시 up   : {lps['up'][:10]}")
    print(f"  예시 down : {lps['down'][:10]}")

    print("\n[공통 유전자 (방향 무관)]")
    common_any = find_common_deg(direction="any")
    print(f"  aging DEG : {common_any['n_aging_deg']}")
    print(f"  lps   DEG : {common_any['n_lps_deg']}")
    print(f"  공통       : {common_any['n_common']}")
    print(f"  랜덤 기대  : {common_any['expected_random']}")
    print(f"  Enrichment: {common_any['enrichment']}x")
    print(f"  Hypergeom p-value: {common_any['p_value']:.2e}")
    print(f"  공통 유전자: {common_any['common_genes'][:20]}")

    print("\n[공통 up (양쪽 다 상승)]")
    common_up = find_common_deg(direction="up")
    print(f"  aging up ∩ lps up: {common_up['n_common']}개")
    print(f"  유전자: {common_up['common_genes']}")

    print("\n[공통 down (양쪽 다 감소)]")
    common_down = find_common_deg(direction="down")
    print(f"  aging down ∩ lps down: {common_down['n_common']}개")
    print(f"  유전자: {common_down['common_genes']}")

    # CSV 로 저장
    df = pd.DataFrame({
        "gene": common_any["common_genes"],
    })
    out_csv = OUTPUTS_DIR / "common_deg_aging_lps.csv"
    df.to_csv(out_csv, index=False)
    print(f"\n[저장] {out_csv}")


# ---------------------------------------------------------------------------
# (3)-2: 공통 DEG spatial 시각화
# ---------------------------------------------------------------------------
def plot_common_deg_spatial(
    genes=None,
    slice_id=0,
    conditions=None,
    point_size=1.5,
    percentile_clip=99.0,
    output_name=None,
):
    """
    공통 DEG 유전자들을 spatial map 위에 발현량으로 시각화.
    Grid: rows=조건, cols=유전자.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    if genes is None:
        genes = ["Gfap", "Apoe", "Cdkn2a", "C4b", "Nfkbia", "Cd47"]
    if conditions is None:
        conditions = [("merfish_control", "90wk"), ("merfish_lps", "24wk")]

    n_rows = len(conditions)
    n_cols = len(genes)
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(n_cols * 3.5, n_rows * 3.5),
                              squeeze=False)

    for i, (dataset, age) in enumerate(conditions):
        adata_raw = subset(dataset, age=age, slices=[slice_id])
        if adata_raw.n_obs == 0:
            for j in range(n_cols):
                axes[i, j].set_title(f"{dataset} {age}: no data")
                axes[i, j].axis("off")
            continue

        donor_counts = adata_raw.obs["donor_id"].value_counts()
        top_donor = str(donor_counts.index[0])
        adata_raw = adata_raw[
            adata_raw.obs["donor_id"].astype(str) == top_donor
        ].copy()

        adata = _prepare_for_deg(adata_raw)

        x = adata.obs["center_x"].astype(float).values
        y = adata.obs["center_y"].astype(float).values

        var_symbols = adata.var["feature_name"].astype(str).values
        symbol_to_idx = {s: i for i, s in enumerate(var_symbols)}

        for j, gene in enumerate(genes):
            ax = axes[i, j]

            if gene not in symbol_to_idx:
                ax.set_title(f"{gene}: not in panel")
                ax.axis("off")
                continue

            gene_idx = symbol_to_idx[gene]
            expr = adata.X[:, gene_idx]
            if hasattr(expr, "toarray"):
                expr = expr.toarray().flatten()
            else:
                expr = np.asarray(expr).flatten()

            vmax = float(np.percentile(expr, percentile_clip)) \
                if expr.max() > 0 else 1.0

            sc = ax.scatter(x, y, c=expr, s=point_size, cmap="viridis",
                            vmin=0, vmax=vmax, linewidths=0)
            ax.set_aspect("equal")
            ax.invert_yaxis()
            ax.set_xticks([])
            ax.set_yticks([])

            if i == 0:
                ax.set_title(gene, fontsize=13, fontweight="bold")
            if j == 0:
                short_donor = top_donor.replace("MsBrainAgingSpatialDonor_", "d")
                ax.set_ylabel(f"{dataset} {age} ({short_donor})", fontsize=10)

            plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)

    fig.suptitle(
        f"Common DEG spatial expression - slice {slice_id} "
        f"(Aging: control 90wk vs LPS: lps 24wk)",
        fontsize=13, y=1.01
    )
    plt.tight_layout()

    if output_name is None:
        output_name = f"common_deg_spatial_slice{slice_id}.png"
    out_path = PROJECT_ROOT / "figures" / output_name
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return str(out_path.relative_to(PROJECT_ROOT))
