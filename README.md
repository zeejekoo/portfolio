# BioNexus Portfolio — 3-Part Bioinformatics + Agent Demo

바이오넥서스 지원용 포트폴리오. 시험 재현 2건 + Agentic 데모 1건.

## 프로젝트 구성

### [task1_bulk_rnaseq/](./task1_bulk_rnaseq/)
CRISPLD2 논문(Himes 2014) bulk RNA-seq 파이프라인 재현.
cutadapt → STAR → featureCounts → DESeq2 → Enrichr.

### [task2_spatial_aging/](./task2_spatial_aging/)
Allen et al. Cell 2023 노화 뇌 MERFISH 재분석.
scanpy · leiden 클러스터링 · cell-type 공간 지도 · aging vs LPS 비교.

### [task3_agent_demo/](./task3_agent_demo/)
위 두 프로젝트 결과를 자연어로 질의하는 agent + Streamlit UI.
Claude API tool-calling으로 CoScientist류 축소판 데모.

## 왜 이런 구성인가

- **Task 1 (bulk RNA-seq)**: 전통적 오믹스 파이프라인 감각 증명
- **Task 2 (MERFISH spatial)**: 최신 spatial · 멀티모달 다루는 감각 증명
- **Task 3 (agent 통합)**: 위 두 개를 재사용 가능한 tool로 묶어서
  agentic 인터페이스 붙임 — 회사가 지금 집중하는 Agentic AI 흐름을
  직접 재현한 미니 CoScientist

Git 커밋 히스토리 = 데이터 세팅 → 각 task 분석 → tool화 → agent →
UI 순서로 쌓아 성장 궤적을 시각화.
