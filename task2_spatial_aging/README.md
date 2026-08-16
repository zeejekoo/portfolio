# Task 2 — Mouse Spatial Brain Atlas 공간전사체 재분석

**논문**: [Allen et al., *Cell* 2023](https://www.cell.com/cell/fulltext/S0092-8674(22)01523-9)  
"Molecular and spatial signatures of mouse brain aging at single-cell resolution"  
(Zhuang Lab, Harvard · DOI 10.1016/j.cell.2022.12.010)

**프로젝트 목표**  
논문 저자의 결론 **"자연 노화 뇌 = 만성 저강도 신경염증 상태로 진행"** 을 원본 데이터로부터 재분석·정량 검증한다. Unsupervised clustering (Leiden) 이 supervised cell type 라벨보다 세밀한 서브타입을 잡는지, 노화와 급성염증(LPS) 사이 공통 유전자 프로그램이 통계적으로 유의미한지, spatial 정보 위에서 동일 anatomical 위치에 재현되는지를 확인한다.

## 데이터

[CELLxGENE collection 31937775](https://cellxgene.cziscience.com/collections/31937775-0602-4e52-a799-b6acdd2bac2e) 원본 그대로 사용 (서브셋 아님):

- `merfish_control.h5ad` — MERFISH 무처치 (378,918 cells × 374 genes, 12 donors)
- `merfish_lps.h5ad` — MERFISH LPS 처리 (345,934 cells × 374 genes, 8 donors)
- `snrnaseq.h5ad` — snRNA-seq 전사체 전체 (79,667 cells × 20,929 genes)

**환경**: conda env `demo` (Python 3.11 · scanpy 1.11.5 · anndata · leidenalg · scikit-learn)

## 분석 항목 (진행 체크)

### (0) 사전 검증 — 저자 워크플로우 확인 (자체 추가)
- [x] MERFISH 프로브 374개가 snRNA-seq HVG 에서 선택됐는지 hypergeom test  
      → **Enrichment 4x, p=6e-85, verdict=conclusive** (`src/gene_overlap.py`)

### (1) 공간 데이터 전처리
- [x] **(1)-1** 두 h5ad 데이터에 `donor` 정보로 개별 슬라이드 식별, 각 샘플별 Cell-type Map (12장)
- [x] **(1)-2** Leiden 클러스터링 (control 23 · LPS 24 clusters)

### (2) 공간 분석 (`merfish_control`)
- [x] **(2)-1** Leiden vs 해부학적 구조 대조 (cortex 6층 · corpus callosum · striatum · hippocampus 재현)
- [x] **(2)-2** cell_type vs Leiden 일치도 (ARI=0.559, NMI=0.814, homogeneity=0.973)
- [x] **(2)-3** Young vs Old microglia 공간분포·비율 (U-shape 패턴 발견)

### (3) 데이터간 비교
- [x] **(3)-1** 노화(Old) vs LPS 공통 발현 변화 유전자 선별 (139 공통, p=2.45e-05)
- [x] **(3)-2** 선별 유전자 Cell-type map 시각화 (Gfap, Apoe, Cdkn2a, C4b, Nfkbia, Cd47)

**전체 결과 리포트**: [REPORT.md](./REPORT.md)

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
주요 dtype 처리 이슈:
- `slice` 컬럼이 카테고리(문자열) 로 저장돼 있어 `subset(slices=[0])` 호출 시 내부 str 캐스팅 필요 → `_isin_str` 헬퍼로 처리
- `center_x/y` 도 `merfish_control` 에서만 카테고리(문자열) 로 저장돼 있음 (LPS 는 float64) → `load()` 에서 자동 float 변환

## 코드 구성

task2_spatial_aging/
├── scripts/
│ ├── download_data.sh # CELLxGENE 3개 h5ad 재현 가능한 다운로드
│ └── build_task2_notebook.py # (선택) ipynb 빌더
├── src/
│ ├── data_loader.py # load / subset / summary + dtype 픽스
│ ├── tools.py # summary / composition / spatial map (color_by='cell_type'|'leiden')
│ ├── gene_overlap.py # (0): MERFISH-HVG hypergeom test
│ ├── clustering.py # (1)-2, (2)-2: Leiden + ARI/NMI + 혼동행렬
│ ├── aging_analysis.py # (2)-3: 나이별 cell type 비율·통계
│ └── comparison.py # (3)-1, (3)-2: DEG + 공통 유전자 + spatial 시각화
├── figures/ # PNG 산출물 (gitignored)
├── outputs/ # 캐시된 h5ad / DEG CSV (gitignored)
└── *.h5ad # 원본 데이터 (gitignored, ~3.1GB)

## 데이터 검증 노트

- **Donor ID `_12` 두 파일 중복**: `MsBrainAgingSpatialDonor_12` 가 control 과 LPS 모두에 존재. 검증 결과 셀 수 (control 33,241 vs LPS 46,179), slice 별 값 크게 다름 → **다른 마우스, ID 재사용**. 처리 방침: 별개 개체로 취급 (paired 분석 불가).

## 재현 방법

```bash
# 환경
export PATH="/BiO/kbioman/kbiomanuser6/miniconda/envs/demo/bin:$PATH"

# 데이터 (없으면 다운로드, 있으면 skip)
bash scripts/download_data.sh

# 저자 워크플로우 검증 (수 초)
python -m src.gene_overlap

# 개별 모듈 자기검증
python -m src.data_loader
python -m src.tools

# Leiden 클러스터링 (약 5~10분, 결과 outputs/ 에 캐시)
python -m src.clustering

# 나이별 microglia 비교
python -m src.aging_analysis

# 공통 DEG (약 2~5분)
python -m src.comparison
```

## 진행 로그

Git 커밋 접두어로 각 단계 구분:
`chore:` 세팅 · `data:` 데이터 · `feat:` 새 기능 · `fix:` 픽스 · `docs:` 문서 · `refactor:` 개편 · `exp:` 노트북
