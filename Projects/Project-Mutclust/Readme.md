# 🧬 MutClust: A Clustering-Based Algorithm for Detecting SARS-CoV-2 Mutation Hotspots

> **Genome-based mutation hotspot detection and severity prediction in COVID-19**

---

## 📅 Project Duration  
**March 2023 – February 2025**  
(Master’s research project, 1st-author publication)

## 🧰 Technology Stack  
- **Languages & Libraries**: Python, Pandas, NumPy, Scikit-learn  
- **Genomic Analysis**: MAFFT, DBSCAN, Shannon Entropy  
- **Immunological Analysis**: HLA Binding Affinity, edgeR, Network Propagation  
- **Databases**: GISAID, STRING, BioGRID

---

## 🎯 Project Overview  

This project aims to develop a **SARS-CoV-2 mutation clustering algorithm (MutClust)** to detect mutation hotspots associated with **COVID-19 severity**, and to investigate their functional relevance by integrating **multi-omics patient data**.

---

## 👨‍💻 Key Contributions  

1. Designed and implemented the **MutClust** algorithm to cluster viral mutations based on entropy and frequency.  
2. Identified **severity-associated mutation hotspots** using patient-derived viral sequences and statistical tests.  
3. Conducted multi-omics integration to explore **immune mechanisms** associated with high-risk mutation signatures.

---

## 🏆 Major Achievements  

### ✅ 1. Development of the MutClust Algorithm  

- **Dataset**: GISAID SARS-CoV-2 genome sequences (N = 224,318, collected between 2020–2022)
- **Key Methods**:
  - **H-score**: Mutation importance score based on Shannon entropy  
  - **MutClust**: A custom DBSCAN-based 1D clustering algorithm using sliding windows
- **Outcomes**:
  - Identified **477 mutation hotspots**
  - Successfully captured critical hotspots (e.g., `c315`, `c442`) missed by frequency-based methods
  - Demonstrated that **H-score** outperformed existing entropy and frequency-based approaches

---

### ✅ 2. Identification of 28 Severity-Associated Hotspots & Immune Mechanism Integration  

- **Clinical Multi-Omics Cohort**: N = 387  
  - Hospitals: *Chungnam National University Hospital*, *Seoul Medical Center*, *Samsung Medical Center*  
  - Data Types: `scRNA-seq`, `patient-derived viral sequences`, `cytokine panel`
  
- **Techniques Used**:  
  `Scikit-learn`, `ANOVA`, `t-test`, `SelectKBest`, `edgeR (DEG analysis)`, `Network Propagation`

- **Key Findings**:
  - Discovered **28 mutation hotspots** associated with clinical severity
  - Defined mutation signatures such as **c315** and **c442** specific to severe patient groups
  - **Cytokine level analysis** showed increased IFNG, TNF, IL-1B, IL-18 in patients with critical mutations
  - **Network propagation** revealed innate immune pathway hubs (e.g., TLR, RIG-I) linked to hotspot-associated viral proteins
  - **DEG analysis from scRNA-seq** showed elevated inflammatory and cytotoxic gene expression in **NK cells, monocytes, and dendritic cells** of hotspot-positive patients
  - Collectively identified **NK cell-mediated immune responses** as key drivers of COVID-19 severity

---

## 📌 Visual Summary  
> (Insert key figures or visualizations: clustering pipeline, H-score heatmap, cytokine profile plots, etc.)

---

## 📄 Related Publication  
> Manuscript under review (1st author) — DOI/link to be added soon.

---

## 🙋‍♀️ Contact  
- Researcher: **Sohyun Yoon**  
- Email: [sohyun.yoon@example.com](mailto:sohyun.yoon@example.com)  
- GitHub: [github.com/sohyun-yoon](https://github.com/sohyun-yoon)

---

