# Aging Brain Spatial Agent (BioNexus Demo)

**목적**: Zhuang Lab의 mouse brain aging MERFISH + snRNA-seq 데이터를 활용해
"자연어 질의 → agent가 tool 호출 → spatial map/셀 조성/유전자 비교 반환" 파이프라인 구축.

## 원본 논문
- Allen et al., *Cell* 2023 — "Molecular and spatial signatures of mouse brain aging at single-cell resolution"
- Data: [CELLxGENE collection](https://cellxgene.cziscience.com/collections/31937775-0602-4e52-a799-b6acdd2bac2e)
- Code (원저자): https://github.com/ZhuangLab/SpatialBrainAgingCell22

## 데이터
- `merfish_control.h5ad` — MERFISH 무처치 (378,918 cells × 374 genes)
- `merfish_lps.h5ad` — MERFISH LPS 처리
- `snrnaseq.h5ad` — snRNA-seq (전사체 전체)

## 환경
Python 3.11 (conda env: `demo`) — scanpy 1.11.5, anndata, leidenalg 등.

## 진행 로그
Git 커밋 히스토리로 관리 — 각 커밋이 한 단계(데이터 준비 → EDA → tool 스켈레톤 → agent → UI).
