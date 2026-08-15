# Task 3 — Nexus Agent Demo (통합)

Task 1(bulk RNA-seq)과 Task 2(spatial transcriptomics) 결과를
자연어로 질의할 수 있는 agentic 인터페이스.

## 구성
- `src/tools.py`: task1·task2 결과를 조회하는 tool 함수 모음
- `src/agent.py`: Claude API tool-calling 루프
- `app.py`: Streamlit UI (질의 창 + 결과 렌더)

## 데모 시나리오 예시
- "old 마우스 slice 0에서 microglia 어디에 몰려있어?"
- "DEX 처리로 발현이 오른 유전자 상위 10개랑 pathway 요약해줘"
- "노화 뇌와 LPS 뇌에서 공통으로 변한 유전자 리스트 만들어줘"

_(Task 1, 2 완료 후 착수)_
