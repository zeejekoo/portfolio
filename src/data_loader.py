"""
노화 뇌 MERFISH + snRNA-seq 데모의 데이터 로딩 유틸리티.

프로젝트 안의 모든 h5ad 접근은 이 모듈을 통해 이루어진다.
경로, 캐싱, 서브셋 규칙을 한 곳에서 관리하기 위함이며,
tool 함수/노트북/Streamlit UI는 scanpy.read_h5ad 를 직접 호출하지
말고 반드시 여기의 load()/subset()/summary() 를 사용한다.

데이터 스키마 (sanity check로 확인함):
    obs:
        - age         : {'4wk','24wk','90wk'} — 카테고리(문자열)
        - donor_id    : 12마리 마우스, 나이 그룹당 약 4마리
        - slice       : {'0','1','2'} — 카테고리(문자열!)
        - cell_type   : 13개 세포 타입 (neuron, oligodendrocyte, ...)
        - center_x, center_y : 세포 중심 좌표 (µm)
    obsm:
        - 'spatial_coords' / 'X_spatial_coords' : 2D 공간 좌표
        - 'X_umap' : 미리 계산된 UMAP
        - 'X_pca'  : 미리 계산된 PCA

dtype 처리 노트:
    CELLxGENE 는 obs 컬럼을 pandas Categorical (문자열 카테고리)로 저장한다.
    'slice' 는 값이 '0','1','2' 처럼 정수처럼 보이지만 실제로는 문자열이므로,
    subset() 은 필터 값과 컬럼을 모두 str로 캐스팅한 뒤 isin() 을 수행한다.
    이 덕분에 호출자는 slices=[0] 이나 slices=['0'] 어느 쪽을 써도 동작한다.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Literal, Optional, Sequence

import anndata as ad
import scanpy as sc


# ---------------------------------------------------------------------------
# 파일 경로 — 데이터를 옮길 때는 여기만 수정하면 된다
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_FILES = {
    "merfish_control": PROJECT_ROOT / "merfish_control.h5ad",
    "merfish_lps":     PROJECT_ROOT / "merfish_lps.h5ad",
    "snrnaseq":        PROJECT_ROOT / "snrnaseq.h5ad",
}

Dataset = Literal["merfish_control", "merfish_lps", "snrnaseq"]
Age = Literal["4wk", "24wk", "90wk"]


# ---------------------------------------------------------------------------
# 로딩 — 각 h5ad 파일을 프로세스 수명 동안 캐시
# ---------------------------------------------------------------------------
@lru_cache(maxsize=3)
def load(dataset: Dataset) -> ad.AnnData:
    """
    지정한 h5ad 파일을 읽고 반환. 프로세스 안에서 캐시된다.

    첫 호출은 파일 크기에 따라 5~30초 소요. 이후 호출은 즉시 반환.
    반환된 AnnData 를 절대 in-place 로 수정하지 말 것.
    (필요하면 subset() 로 복사본을 받아서 수정)
    """
    path = DATA_FILES[dataset]
    if not path.exists():
        raise FileNotFoundError(
            f"{path} 파일이 없습니다. 먼저 scripts/download_data.sh 를 실행하세요."
        )
    return sc.read_h5ad(path)


# ---------------------------------------------------------------------------
# 서브셋 — 자주 걸리는 4개 축의 필터를 하나로
# ---------------------------------------------------------------------------
def _isin_str(series, values) -> "pd.Series[bool]":
    """
    dtype-무관 isin: 컬럼과 값을 모두 str 로 캐스팅한 뒤 비교.
    Categorical / object / 숫자 컬럼 어디에서든 안전.
    """
    return series.astype(str).isin([str(v) for v in values])


def subset(
    dataset: Dataset,
    age: Optional[Age | Sequence[Age]] = None,
    cell_types: Optional[Sequence[str]] = None,
    donor_ids: Optional[Sequence[str]] = None,
    slices: Optional[Sequence[int | str]] = None,
) -> ad.AnnData:
    """
    주어진 나이/세포타입/donor/슬라이스 조건으로 잘라낸 AnnData 복사본을 반환.
    None 인 인자는 해당 축에 필터를 걸지 않는다는 뜻.

    반환값은 항상 사본이므로 자유롭게 수정해도 원본에는 영향 없음.
    """
    adata = load(dataset)
    mask = None

    def _combine(m, new):
        return new if m is None else (m & new)

    if age is not None:
        ages = [age] if isinstance(age, str) else list(age)
        mask = _combine(mask, _isin_str(adata.obs["age"], ages).values)
    if cell_types is not None:
        mask = _combine(mask, _isin_str(adata.obs["cell_type"], cell_types).values)
    if donor_ids is not None:
        mask = _combine(mask, _isin_str(adata.obs["donor_id"], donor_ids).values)
    if slices is not None:
        mask = _combine(mask, _isin_str(adata.obs["slice"], slices).values)

    if mask is None:
        return adata.copy()
    return adata[mask].copy()


# ---------------------------------------------------------------------------
# 요약 — 셀 수·나이·donor·세포타입 상위 목록을 dict 로
# ---------------------------------------------------------------------------
def summary(dataset: Dataset) -> dict:
    """agent 컨텍스트 / 상태바 / 로그용 압축 요약 dict."""
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


# ---------------------------------------------------------------------------
# 자기 검증 — `python -m src.data_loader` 로 실행
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import json
    s = summary("merfish_control")
    s["cell_types"] = dict(list(s["cell_types"].items())[:5]) | {"...": "..."}
    print(json.dumps(s, indent=2, default=str, ensure_ascii=False))

    # dtype 버그 회귀 방지 체크
    n = subset("merfish_control", age="90wk", slices=[0]).n_obs
    print(f"\n회귀 체크: age=90wk, slices=[0] → {n} cells (0 이상이어야 정상)")
