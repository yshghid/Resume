## Mutclust: 돌연변이 핫스팟 탐지 알고리즘 

- 수행기간: 2023.03 – 2025.02 (석사 연구, 제1저자 논문)

- 사용 기술: Python, Pandas, NumPy, Scikit-learn, MAFFT, edgeR, STRING, Biogrid, DBSCAN, Shannon Entropy, HLA Binding Affinity Analysis, Network Propagation

- 프로젝트 정의 및 목표: COVID-19 중증도 예측을 위한 **SARS-CoV-2 유전체 분석 알고리즘** 개발

- 주요 역할
  1) 돌연변이 클러스터링 알고리즘(MutClust) 고안 및 구현
  2) 중증도와 연관된 돌연변이 hotspot 탐지 및 Multi-omics 환자군 데이터 연동 분석 수행
  3) 중증 연관 hotspot 보유군의 면역학적 기전 규명

- 핵심 성과
  1) Mutation Hotspot 탐지를 위한 MutClust 알고리즘 개발
     - 데이터: GISAID SARS-CoV-2 유전체 시퀀스 (N=224,318/2020–2022)
     - 사용 스택 및 기술: Python, NumPy, Pandas, H-score(Shannon Entropy 기반 돌연변이 중요도), MutClust(In-house 유전체 분석 알고리즘)
     - 성과
       1) 빈도 + 다양성 기반 1D 클러스터링으로 477개 mutation hotspot 탐지
       2) 기존 빈도 기반 방식이 놓친 중증 관련 핫스팟(c315, c442) 포착 가능
       3) H-중요도가 기존 nucleotide entropy 및 빈도 기반 접근법보다 높은 설명력을 보임

  2) 중증 관련 돌연변이 Hotspot 28개 도출 및 Multi-omics 통합 기전 구축
     - 데이터
       - COVID-19 multi-omics 코호트(N=387/Chungnam National University Hospital, Seoul Medical Center, Samsung Medical Center)
       - scRNA-seq, patient driven SARS-CoV-2 sequence, Cytokine
     - 사용 스택 및 기술: Scikit-learn (SelectKBest, Clustermap), ANOVA, t-test, edgeR(DE analysis)
     - 성과
       1) 중증 관련 돌연변이 Hotspot 28개 도출
       2) Severity 연관 hotspot 보유군의 mutation signature(c315, c442) 도출
       3) 사이토카인 수준 비교, network propagation 분석, DEG 분석 결과
          - 중증 hotspot 보유 환자군에서 염증성 사이토카인(IFNG, TNF 등)의 증가와 NK cell에서의 염증 및 세포독성 관련 유전자 발현 증가 확인
          - COVID-19 중증 예후의 NK 세포 연관 면역학적 기전 규명


