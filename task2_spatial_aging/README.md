# Task 2 — Mouse Spatial Brain Atlas (공간전사체 재분석)

**논문**: [Allen et al., *Cell* 2023](https://www.cell.com/cell/fulltext/S0092-8674(22)01523-9)  
"Molecular and spatial signatures of mouse brain aging at single-cell resolution"  
(Zhuang Lab, Harvard · DOI 10.1016/j.cell.2022.12.010)

**데이터**: [CELLxGENE 컬렉션](https://cellxgene.cziscience.com/collections/31937775-0602-4e52-a799-b6acdd2bac2e)
- `merfish_control.h5ad` — MERFISH 무처치 (378,918 cells × 374 genes)
- `merfish_lps.h5ad`     — MERFISH LPS 처리 (345,934 cells × 374 genes)
- `snrnaseq.h5ad`        — snRNA-seq 전사체 전체 (79,667 cells × 20,929 genes)

**환경**: conda env `demo` (Python 3.11 + scanpy 1.11.5)

---

## 시험 요구사항 (진행 체크리스트)

### (0) 사전 검증 — 저자 워크플로우 확인 (자체 추가)
- [x] MERFISH 프로브 374개가 snRNA-seq HVG 에서 선택됐는지 hypergeom test  
      → **Enrichment 4x, p=6e-85, verdict=conclusive** (`src/gene_overlap.py`)

### (1) 공간 데이터 전처리
- [ ] **(1)-1** 두 h5ad 데이터에 `donor` 정보로 개별 슬라이드 식별,  
      각 샘플별 Cell-type Map 그리기
- [ ] **(1)-2** Leiden 클러스터링 수행

### (2) 공간 분석 (`merfish_control` 데이터만)
- [ ] **(2)-1** Leiden 클러스터가 뇌 조직의 해부학적 구조  
      (cortex 층 등) 를 얼마나 반영하는지 시각적 대조
- [ ] **(2)-2** 기존 `cell_type` vs Leiden 일치도 평가 (ARI/NMI)  
      + 불일치 클러스터 원인 추정
- [ ] **(2)-3** `age` 활용, Young vs Old 마우스 간 특정 세포 유형  
      (예: Microglia) 의 공간분포·클러스터 비율 변화 정량 비교

### (3) 데이터간 비교
- [ ] **(3)-1** 노화(Old) vs LPS 두 샘플 공통 발현 변화 유전자 선별
- [ ] **(3)-2** 선별 유전자군을 Cell-type map 위에 시각화 비교

---

## 데이터 스키마 (확정)
obs 컬럼:
age {'4wk','24wk','90wk'} 카테고리 문자열
donor_id 12마리 마우스 카테고리 문자열
slice {'0','1','2'} 카테고리 문자열 (!)
cell_type 13종 neuron, oligodendrocyte, microglial cell, ...
center_x/y µm 단위 좌표
obsm:
spatial_coords / X_spatial_coords 2D 물리 좌표
X_umap 미리 계산된 UMAP
X_pca 미리 계산된 PCA (50 dim)
adata.raw:
원본 정수 카운트 살아있음 (재현 가능)
`slice` 는 정수처럼 보이지만 카테고리 문자열이므로 `subset(slices=[0])` 호출 시 내부에서 str 캐스팅 필요. `src/data_loader.py` 의 `_isin_str` 헬퍼가 처리.

## 코드 구성
task2_spatial_aging/
├── scripts/download_data.sh # CELLxGENE 3개 h5ad 재현 가능한 다운로드
├── src/
│ ├── data_loader.py # load / subset / summary + dtype 픽스
│ ├── tools.py # summary / composition / spatial map
│ ├── gene_overlap.py # (0) 검증: hypergeom test
│ ├── clustering.py # (예정) (1)-2, (2)-2: Leiden + ARI/NMI
│ └── comparison.py # (예정) (3)-1: aging vs LPS
├── notebooks/
│ └── task2_full_analysis.ipynb # (예정) 시험 제출 형식과 일치
├── figures/ # PNG 산출물 (gitignored)
├── outputs/ # 캐시된 h5ad (gitignored)
└── *.h5ad # 원본 데이터 (gitignored, ~3.1GB)
## 재현 방법

```bash
export PATH="/BiO/kbioman/kbiomanuser6/miniconda/envs/demo/bin:$PATH"
bash scripts/download_data.sh
python -m src.gene_overlap
python -m src.data_loader
python -m src.tools
```

## 진행 로그

Git 커밋 접두어로 각 단계 구분:
`chore:` 세팅 · `data:` 데이터 · `feat:` 새 기능 · `fix:` 픽스 · `docs:` 문서 · `refactor:` 개편 · `exp:` 노트북
