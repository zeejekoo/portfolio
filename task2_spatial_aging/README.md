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


## 재현 방법

```bash
# 환경
export PATH="/BiO/kbioman/kbiomanuser6/miniconda/envs/demo/bin:$PATH"

# 데이터 (없으면 다운로드, 있으면 skip)
bash scripts/download_data.sh

# 저자 워크플로우 검증
python -m src.gene_overlap

# 개별 tool 자기검증
python -m src.data_loader
python -m src.tools
```

## 진행 로그

Git 커밋 접두어로 각 단계 구분:
- `chore:` 세팅   `data:` 데이터 스크립트   `feat:` 새 기능
- `fix:` 버그 픽스   `docs:` 문서   `refactor:` 구조 개편
- `exp:` 실험·EDA·노트북
