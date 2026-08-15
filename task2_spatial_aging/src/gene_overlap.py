"""
데이터셋 간 유전자 겹침 검증 (hypergeometric test).

주 사용 시나리오:
    MERFISH panel(374) 이 snRNA-seq HVG(2958) 안에서 선택됐는지 검증.
    저자의 논문 워크플로우(snRNA-seq 로 marker 발견 → MERFISH 프로브 
    디자인) 가 실제 데이터에서 관찰되는지 정량 확인.

세 개 함수:
    get_gene_names()     : 데이터셋의 유전자 이름 집합 (선택: HVG만)
    enrichment_test()    : 범용 hypergeometric enrichment test
    check_probe_selection(): 위 둘을 엮은 편의 래퍼 (판정 verdict 포함)

Task 3 agent tool 로 그대로 노출 예정.
"""

from __future__ import annotations

import anndata as ad
import scanpy as sc
from scipy.stats import hypergeom

from src.data_loader import Dataset, load


def get_gene_names(dataset: Dataset, hvg_only: bool = False) -> set:
    """
    데이터셋의 유전자 이름 집합을 반환.

    hvg_only=True 면 var['highly_variable'] 이 True 인 유전자만.
    MERFISH 같은 targeted panel 은 HVG 컬럼이 없으므로 이 경우 전체 반환.
    """
    adata = load(dataset)
    if hvg_only and "highly_variable" in adata.var.columns:
        mask = adata.var["highly_variable"].astype(bool)
        genes = adata.var.loc[mask, "feature_name"].astype(str)
    else:
        genes = adata.var["feature_name"].astype(str)
    return set(genes)


def enrichment_test(
    sample: set,
    success_pool: set,
    universe: set,
) -> dict:
    """
    Hypergeometric enrichment test — 범용.

    Universe 안에서 success_pool 이 차지하는 비율을 기준선으로,
    sample 안에서 success_pool 원소가 얼마나 과대표현됐는지 계산.

    파라미터:
        sample       : 관찰 집합 (예: MERFISH 374 유전자)
        success_pool : 관심 집합 (예: snRNA HVG 2958 유전자)
        universe     : 전체 배경 (예: snRNA 20926 유전자)

    반환 dict:
        n_universe  : 전체 크기 (M)
        n_success   : 성공 pool 크기 (n)
        n_sample    : 표본 크기 (N, universe 로 제한 후)
        observed    : 실제 겹침 (k)
        expected    : 랜덤 기대값 (N × n/M)
        enrichment  : effect size = observed / expected
        p_value     : P(X >= observed) — hypergeometric survival function
    """
    # universe 로 제한 (외부 유전자 있으면 제거)
    sample = sample & universe
    success_pool = success_pool & universe

    M = len(universe)
    n = len(success_pool)
    N = len(sample)
    k = len(sample & success_pool)

    expected = N * n / M if M > 0 else 0.0
    enrichment = k / expected if expected > 0 else float("inf")
    p_value = float(hypergeom.sf(k - 1, M, n, N))

    return {
        "n_universe": M,
        "n_success": n,
        "n_sample": N,
        "observed": k,
        "expected": round(expected, 2),
        "enrichment": round(enrichment, 2),
        "p_value": p_value,
    }


def check_probe_selection(
    panel_dataset: Dataset = "merfish_control",
    pool_dataset: Dataset = "snrnaseq",
    use_hvg: bool = True,
) -> dict:
    """
    'panel 이 pool 에서 유의미하게 선택됐는가' 를 판정하는 편의 함수.

    기본값: MERFISH(374 targeted panel) 이 snRNA-seq HVG(2958) 에서
    선택됐는가? → hypergeom test + verdict 문자열 반환.

    반환값 추가 필드:
        panel_dataset / pool_dataset / pool_type : 인자 기록
        verdict : "random-like" / "moderate" / "strong" / "conclusive"
    """
    panel_genes = get_gene_names(panel_dataset)              # sample
    pool_all = get_gene_names(pool_dataset, hvg_only=False)  # universe
    pool_success = get_gene_names(pool_dataset, hvg_only=use_hvg)  # success

    result = enrichment_test(panel_genes, pool_success, pool_all)

    # p-value 기반 자동 판정
    p = result["p_value"]
    if p > 0.05:
        verdict = "random-like"
    elif p > 1e-10:
        verdict = "moderate"
    elif p > 1e-50:
        verdict = "strong"
    else:
        verdict = "conclusive"

    result.update({
        "panel_dataset": panel_dataset,
        "pool_dataset": pool_dataset,
        "pool_type": "HVG" if use_hvg else "all_genes",
        "verdict": verdict,
    })
    return result


if __name__ == "__main__":
    print("=== MERFISH panel vs snRNA-seq HVG 겹침 검증 ===\n")
    result = check_probe_selection(
        panel_dataset="merfish_control",
        pool_dataset="snrnaseq",
        use_hvg=True,
    )

    print(f"Universe (snRNA-seq total)  : {result['n_universe']:>8d}")
    print(f"Success pool (snRNA HVG)    : {result['n_success']:>8d}")
    print(f"Sample (MERFISH panel)      : {result['n_sample']:>8d}")
    print()
    print(f"Observed overlap            : {result['observed']:>8d}")
    print(f"Expected (random)           : {result['expected']:>8.1f}")
    print(f"Enrichment (effect size)    : {result['enrichment']:>8.2f}x")
    print(f"P-value (hypergeometric)    : {result['p_value']:.2e}")
    print(f"Verdict                     : {result['verdict']}")
