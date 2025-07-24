# 🧬 MutClust: SARS-CoV-2 돌연변이 핫스팟 탐지 알고리즘

> **COVID-19 중증도 예측을 위한 유전체 기반 돌연변이 클러스터링 분석**

---

## 📅 수행 기간  
**2023.03 – 2025.02** (석사 연구, 제1저자 논문 기반)

## 🧰 사용 기술 스택  
- Python, Pandas, NumPy, Scikit-learn, MAFFT, edgeR, STRING, Biogrid, DBSCAN, Shannon Entropy, HLA Binding Affinity Analysis, Network Propagation 

---

## 🎯 프로젝트 개요  

본 프로젝트는 **SARS-CoV-2 유전체의 돌연변이 특성**을 기반으로 COVID-19 **중증도 예측에 기여하는 mutation hotspot을 탐지**하고, 이를 **환자 면역학 데이터와 통합 분석**하여 **질병 악화 기전**을 규명하는 것을 목표로 하였습니다.

---

## 👨‍💻 주요 역할  

1. **MutClust 알고리즘 개발**: 빈도 및 엔트로피 기반의 mutation 중요도(H-score)를 반영한 유전체 클러스터링 알고리즘 설계 및 구현
2. **Multi-omics 환자 데이터 연동 분석**: 중증 환자군의 돌연변이 패턴과 면역 기전 간의 연관성 탐색  

---

## 🏆 핵심 성과

### ✅ 1. MutClust 알고리즘 개발

- **데이터**: GISAID SARS-CoV-2 유전체 시퀀스 (N = 224,318 / 2020–2022)
- **기술 요약**:
  - **H-score**: Shannon Entropy 기반 돌연변이 중요도 정의
  - **MutClust**: Sliding-window 기반 DBSCAN 클러스터링 알고리즘
- **성과**:
  - 총 **477개 mutation hotspot** 탐지
  - 기존 빈도 기반 방법이 놓친 중증 관련 hotspot (**c315**, **c442**) 포착 가능
  - **H-score**가 기존 방법보다 높은 설명력 확보

---

### ✅ 2. 중증 연관 Hotspot 28개 도출 및 면역 기전 분석

- **데이터 출처**:
  - **COVID-19 multi-omics cohort (N = 387)**  
    - *Chungnam National University Hospital*  
    - *Seoul Medical Center*  
    - *Samsung Medical Center*
  - 포함 데이터: `scRNA-seq`, `patient-derived viral sequence`, `cytokine panel`

- **분석 스택**:  
  `Scikit-learn`, `ANOVA`, `t-test`, `edgeR`, `SelectKBest`, `Network Propagation`

- **주요 결과**:
  - 중증 관련 **돌연변이 hotspot 28개** 도출
  - 중증 환자에서 **대표 mutation signature**: `c315`, `c442`
  - **사이토카인 수준 비교**: IFNG, TNF 등의 **염증성 cytokine 증가**
  - **Network propagation**: 선천면역 기반 **숙주-바이러스 상호작용 경로** 활성화
  - **DEG 분석**: NK cell, monocyte 등에서 **염증 및 세포독성 유전자 발현 증가**
  - 최종적으로, **중증 악화와 관련된 NK 세포 중심의 면역학적 기전** 규명

---

## 📄 관련 논문 (In Review)  
> 작성 중, 곧 링크 예정

---

## 🙋‍♀️ Contact  
- 연구자: **윤소현 (Sohyun Yoon)**  
- 이메일: [sohyun.yoon@example.com](mailto:sohyun.yoon@example.com)  
- GitHub: [github.com/sohyun-yoon](https://github.com/sohyun-yoon)

---

