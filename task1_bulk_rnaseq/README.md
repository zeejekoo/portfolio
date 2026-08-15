# Task 1 — CRISPLD2 Bulk RNA-seq Pipeline

**논문**: Himes et al., PLOS One 2014 (doi:10.1371/journal.pone.0099625)
**데이터**: GEO GSE52778 
## 파이프라인
1. cutadapt로 5' 12bp trim
2. STAR 인덱스 빌드 + 매핑
3. Mapping rate 논문값(83.36%) 대비 비교
4. featureCounts + DESeq2로 Control vs DEX DEG
5. Enrichr로 pathway enrichment
6. 논문 결론(CRISPLD2·항염증) 지지 여부 판단

**환경**: conda env `rnaseq` (bioconda: cutadapt, star, samtools, subread)
_(작업 예정)_
