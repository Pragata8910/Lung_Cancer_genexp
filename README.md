# Gene Expression Analysis in Non Small Cell Lung Cancer(NSCLC) Using RNAseq Data.
## Overview
A computational analysis of differential gene expression patterns in non-small cell lung cancer under various parameters using RNA sequencing data.
### Brief Project Description
This project analyzes the raw high-throughput RNA sequencing data to investigate gene expression patterns in non-small cell lung cancers.

### The analysis includes:

* Differential Gene Expression Analysis in NSCLC cells according to various parameters:
  * Smoking status
  * Histology
  * Gender
  * Age
* A specific function has been developed to plot the differential expression of key genes:
  * EGFR
  * BRAF etc.
  * And other significant markers
* Pathway Enrichment Analysis
* Visualization of the key findings
* Statistical validation of the expression patterns

## **Respository Structure**

**RNA-Seq Annalysis**

|---**data**

|      |--count_matrix

|      |--raw & processed metadata

|      |--results for GO & KEGG analysis

|      |--significant genes as outcome of DGE Analysis

|

|---**results(Visualisations from analysis)**

|      |--MA_plot.png

|      |--Volcano_plot.png

|      |--Heatmap.pdf

|      |--Other plots to visualise significant and key genes

|

|---**scripts**

|      |--01_dataprocessing.R(preprocessing of raw data)

|      |--02_Differential_expression.R(DeSeq2 analysis)

|      |--03_Visualisationns.R

|      |--04_Pathway_analysis.R

|

|---**README.md**

|---**LICENSE**

|---**.gitignore**


## Dataset Information
- **Source:** [GEO Accession: GSE81089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81089)  
- **Data Type:** High Throughput RNA Seq data
- **Samples:** Non Small Cell Lung cancer patients (tumor & normal)  
- **Format:** Read counts (raw), Metadata  

##  Methodology
1️**Data Preprocessing**  
   - Loaded raw count data & metadata  
   - Filtered & matched samples between count matrix and metadata  
   - Normalized counts using **DESeq2**  

2️ **Differential_expression (DESeq2)**
   - Performed differential expression with **`DESeqDataSetFromMatrix()`**  
   - Design formula: `~ condition`  
   - Identified significantly up/downregulated genes (`padj < 0.2`)  

3️ **Gene Annotation**
   - Converted **Ensembl IDs → Gene Symbols**  
   - Mapped genes to biological functions  

4️ **Pathway Enrichment Analysis**
   - **GO (Gene Ontology) analysis**  
   - **KEGG Pathway analysis**  

5️ **Visualization**
   - **MA Plot**: Log2FoldChange vs mean expression  
   - **Volcano Plot**: Significant genes based on log2FC & padj  
   - **Heatmap**: Top differentially expressed genes  
   - **Gene Expression Boxplots & Barplots**: Specific genes (EGFR, TP53, ALK) across conditions  

## Key Findings
- One of the most important findings include the significant upregulation of **MMP9** in Squamous cell carcinoma compared to Adenocarinoma which explains the relative difference in invasiveness between the two
- **EGFR, TP53, KRAS, MET, PIK3CA, ALK, ROS1** showed significant expression differences between tumor and normal samples.
- Most of these gens were consistently upregulated in all kinds of tumors subtypes with some variation
- Pathway analysis highlighted **AMPK signalling**,**MAPK signaling**, **p53 signaling**, and **immune response pathways** as critical in NSCLC.  
- Smoking, histology significantly influenced gene expression, altering the expression of multiple oncogenes.  

## How to Reproduce This Analysis
###   **Prerequisites**
- **R (≥4.0)**
- **R Packages:**
  ```r
  install.packages(c("DESeq2", "tidyverse", "pheatmap", "ggplot2", "clusterProfiler", "org.Hs.eg.db")'''
1. **Clone the repository**

 Run the following commands in your terminal stepwise-
 
    git clone https://github.com/Pragata8910/Lung_Cancer_genexp.git
    
    cd RNA_Seq_Analysis

2.**Run the scripts in order**

    source("scripts/data_processing.R)
    
    source("Differential_expression.R")
    
    source("Pathway_analysis.R")
    
    source("Visualisations.R")

3. Check for the plots/result folder

## **License**

This project is licensed under the [MIT License](LICENSE)

### Feel free to fork this repository and open a pull request if you want to contribute

## **Contact**

For any queries, reach out via 
* Email: [Pragata Ghosh](mailto:pragata2004@gmail.com)
* GitHub: [@Pragata8910](https://github.com/Pragata8910)







     
