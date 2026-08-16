"""
Leiden 클러스터링 실행/평가.

(1)-2, (2)-2 요구사항 충족.

함수:
    run_leiden()              — Leiden 실행 + 파일 캐시  ((1)-2)
    evaluate_clustering()     — ARI/NMI/혼동행렬 계산     ((2)-2)
    plot_confusion_heatmap()  — 혼동행렬 히트맵 PNG      ((2)-2)

Task 3 agent 가 이 함수들을 tool 로 노출 예정.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    homogeneity_score,
    completeness_score,
    v_measure_score,
)

from src.data_loader import Dataset, PROJECT_ROOT, load


OUTPUTS_DIR = PROJECT_ROOT / "outputs"
OUTPUTS_DIR.mkdir(exist_ok=True)

FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


def _cache_path(dataset: str, resolution: float) -> Path:
    return OUTPUTS_DIR / f"leiden_{dataset}_r{resolution:g}.h5ad"


# ---------------------------------------------------------------------------
# 1. Leiden 실행 ((1)-2)
# ---------------------------------------------------------------------------
def run_leiden(
    dataset: Dataset,
    resolution: float = 1.0,
    use_rep: str = "X_pca",
    n_neighbors: int = 15,
    force_recompute: bool = False,
) -> ad.AnnData:
    cache = _cache_path(dataset, resolution)
    if cache.exists() and not force_recompute:
        print(f"[cache 히트] {cache.name}")
        return sc.read_h5ad(cache)

    print(f"[계산 시작] {dataset}, resolution={resolution}")
    adata = load(dataset).copy()
    if use_rep not in adata.obsm:
        raise KeyError(f"obsm 에 '{use_rep}' 없음")

    print(f"  [1/3] 이웃 그래프")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)

    print(f"  [2/3] Leiden")
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph",
                 n_iterations=2, directed=False)

    print(f"  [3/3] 캐시 저장: {cache.name}")
    adata.write_h5ad(cache)
    print(f"[완료] {dataset}: {adata.n_obs:,} cells, "
          f"{adata.obs['leiden'].nunique()} clusters")
    return adata


# ---------------------------------------------------------------------------
# 2. 클러스터링 평가 ((2)-2)
# ---------------------------------------------------------------------------
def evaluate_clustering(
    dataset: Dataset,
    resolution: float = 1.0,
    cluster_key: str = "leiden",
    label_key: str = "cell_type",
) -> dict:
    """
    Leiden 클러스터 vs supervised label 정량 비교.

    지표 4종:
        ARI          : Adjusted Rand Index (전체 일치도, -1~1, 1=완전 일치)
        NMI          : Normalized Mutual Information (0~1, 정보량 기반)
        homogeneity  : 각 클러스터 내부 순도 (하나의 cell_type 으로만 이뤄졌나)
        completeness : 각 cell_type 통합도 (하나의 클러스터로 모였나)
        v_measure    : homogeneity·completeness 의 조화평균

    반환:
        위 지표 + 혼동행렬 3종 (raw counts, row-norm, col-norm) 담긴 dict
    """
    cache = _cache_path(dataset, resolution)
    if not cache.exists():
        raise FileNotFoundError(
            f"Leiden 캐시 없음: {cache}. 먼저 run_leiden() 실행."
        )

    adata = sc.read_h5ad(cache)
    labels = adata.obs[label_key].astype(str).values
    clusters = adata.obs[cluster_key].astype(str).values

    ari = float(adjusted_rand_score(labels, clusters))
    nmi = float(normalized_mutual_info_score(labels, clusters))
    hom = float(homogeneity_score(labels, clusters))
    com = float(completeness_score(labels, clusters))
    vm = float(v_measure_score(labels, clusters))

    conf = pd.crosstab(
        pd.Series(labels, name=label_key),
        pd.Series(clusters, name=cluster_key),
    )
    # 열(leiden) 을 숫자 순으로 정렬
    try:
        col_order = sorted(conf.columns, key=lambda x: int(x))
        conf = conf[col_order]
    except ValueError:
        pass

    conf_row_norm = conf.div(conf.sum(axis=1), axis=0).round(3)
    conf_col_norm = conf.div(conf.sum(axis=0), axis=1).round(3)

    return {
        "dataset": dataset,
        "n_cells": int(adata.n_obs),
        "n_labels": int(conf.shape[0]),
        "n_clusters": int(conf.shape[1]),
        "ari": ari,
        "nmi": nmi,
        "homogeneity": hom,
        "completeness": com,
        "v_measure": vm,
        "confusion_counts": conf,
        "confusion_row_norm": conf_row_norm,
        "confusion_col_norm": conf_col_norm,
    }


# ---------------------------------------------------------------------------
# 3. 혼동행렬 히트맵
# ---------------------------------------------------------------------------
def plot_confusion_heatmap(
    dataset: Dataset,
    resolution: float = 1.0,
    normalize: str = "row",   # "row" | "col" | "none"
    output_name: Optional[str] = None,
) -> str:
    """
    혼동행렬 히트맵을 PNG 로 저장. 반환: 파일 경로.

    normalize="row" (기본, 가장 해석 편함):
        각 cell_type (행) 안에서 어느 leiden 클러스터 (열) 로 흩어졌는지 비율.
        1.0 = 그 라벨 전부가 한 클러스터로 모임 (completeness 완벽)
        0.5 씩 두 클러스터에 = 서브타입으로 분리됐음
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    result = evaluate_clustering(dataset, resolution)

    if normalize == "row":
        matrix = result["confusion_row_norm"]
        title_suffix = "row-normalized (라벨별 클러스터 분포)"
        cmap = "YlOrRd"
    elif normalize == "col":
        matrix = result["confusion_col_norm"]
        title_suffix = "col-normalized (클러스터별 라벨 구성)"
        cmap = "YlOrRd"
    else:
        matrix = result["confusion_counts"]
        title_suffix = "raw counts"
        cmap = "Blues"

    fig, ax = plt.subplots(
        figsize=(matrix.shape[1] * 0.4 + 3, matrix.shape[0] * 0.35 + 2)
    )
    im = ax.imshow(matrix.values, cmap=cmap, aspect="auto")

    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=0, fontsize=8)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=9)

    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("cell_type")
    ax.set_title(
        f"{dataset} — cell_type vs Leiden {title_suffix}\n"
        f"ARI={result['ari']:.3f}  NMI={result['nmi']:.3f}  "
        f"V={result['v_measure']:.3f}  "
        f"(hom={result['homogeneity']:.3f}, com={result['completeness']:.3f})"
    )

    # 셀 값 표시 (유의미한 값만)
    if normalize != "none":
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                v = matrix.iloc[i, j]
                if v > 0.1:
                    color = "white" if v > 0.5 else "black"
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            color=color, fontsize=7)

    plt.colorbar(im, ax=ax, fraction=0.03)

    if output_name is None:
        output_name = f"confusion_{dataset}_r{resolution:g}_{normalize}.png"
    out_path = FIGURES_DIR / output_name
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(out_path.relative_to(PROJECT_ROOT))


# ---------------------------------------------------------------------------
# self-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    for ds in ["merfish_control", "merfish_lps"]:
        print("=" * 60)
        print(f"### {ds} ###")
        res = evaluate_clustering(ds)
        print(f"셀 {res['n_cells']:,}, 라벨 {res['n_labels']}, "
              f"클러스터 {res['n_clusters']}")
        print(f"  ARI          : {res['ari']:.3f}")
        print(f"  NMI          : {res['nmi']:.3f}")
        print(f"  homogeneity  : {res['homogeneity']:.3f}  "
              f"(각 클러스터 순도)")
        print(f"  completeness : {res['completeness']:.3f}  "
              f"(각 라벨 통합도)")
        print(f"  V-measure    : {res['v_measure']:.3f}")

        p = plot_confusion_heatmap(ds, normalize="row")
        print(f"  히트맵 저장: {p}")
        print()
