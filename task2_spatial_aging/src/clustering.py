"""
Leiden 클러스터링 실행과 캐싱.

시험 Task 2 (1)-2 요구사항.

핵심 설계:
    - X_pca (저자가 이미 계산) 를 재활용해서 이웃 그래프 계산
    - Leiden 결과를 outputs/leiden_{dataset}_r{res}.h5ad 로 캐시
    - 캐시 있으면 즉시 로드, 없으면 계산 후 저장

Task 3 의 agent 가 이 함수를 그대로 tool 로 노출 예정.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional

import anndata as ad
import scanpy as sc

from src.data_loader import Dataset, PROJECT_ROOT, load


OUTPUTS_DIR = PROJECT_ROOT / "outputs"
OUTPUTS_DIR.mkdir(exist_ok=True)


def _cache_path(dataset: str, resolution: float) -> Path:
    return OUTPUTS_DIR / f"leiden_{dataset}_r{resolution:g}.h5ad"


def run_leiden(
    dataset: Dataset,
    resolution: float = 1.0,
    use_rep: str = "X_pca",
    n_neighbors: int = 15,
    force_recompute: bool = False,
) -> ad.AnnData:
    """
    지정한 데이터셋에 Leiden 클러스터링 수행.
    obs['leiden'] 이 채워진 AnnData 반환.

    파라미터:
        dataset          : merfish_control / merfish_lps / snrnaseq
        resolution       : Leiden 해상도 (기본 1.0)
        use_rep          : 이웃 계산에 쓸 임베딩 (기본: X_pca)
        n_neighbors      : k-NN 그래프 k (기본 15)
        force_recompute  : True 면 캐시 무시하고 재계산
    """
    cache = _cache_path(dataset, resolution)

    if cache.exists() and not force_recompute:
        print(f"[cache 히트] {cache.name}")
        return sc.read_h5ad(cache)

    print(f"[계산 시작] {dataset}, resolution={resolution}")
    adata = load(dataset).copy()

    # X_pca 존재 확인
    if use_rep not in adata.obsm:
        raise KeyError(
            f"obsm 에 '{use_rep}' 없음. 저자가 배포한 PCA 결과가 필요."
        )

    print(f"  [1/3] 이웃 그래프 계산 (n_neighbors={n_neighbors}, use_rep={use_rep})")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)

    print(f"  [2/3] Leiden 실행 (resolution={resolution}, flavor=igraph)")
    sc.tl.leiden(
        adata,
        resolution=resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    print(f"  [3/3] 캐시 저장: {cache.name}")
    adata.write_h5ad(cache)

    n_clusters = adata.obs["leiden"].nunique()
    print(f"[완료] {dataset}: {adata.n_obs:,} cells, {n_clusters} clusters")
    return adata


if __name__ == "__main__":
    # 시험 (1)-2: 두 MERFISH 데이터셋에 Leiden 수행
    for ds in ["merfish_control", "merfish_lps"]:
        print("=" * 60)
        adata = run_leiden(ds, resolution=1.0)
        print(f"\n클러스터별 셀 수 (상위 10):")
        print(adata.obs["leiden"].value_counts().head(10))
        print()
