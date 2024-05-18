# DRN-AD-EWAS
# Analytical Pipeline for DNA Modification Profiling in the Dorsal Raphe Nucleus

This repository contains the R scripts used for the analytical pipeline described in the manuscript "DNA modification profiling in the dorsal raphe nucleus highlights cell type-specific changes in TNXB in Alzheimer’s disease" (Riemens et al., 2024).
Overview

In this study, we conducted an epigenome-wide association study (EWAS) to profile DNA methylation and hydroxymethylation in the dorsal raphe nucleus (DRN) of post-mortem brain tissue from individuals with Alzheimer's disease (AD). The repository includes all the necessary R code to reproduce the analysis, from data preprocessing to differential modification analysis, and gene ontology enrichment.
Contents

Data Preprocessing: 
Scripts for loading and preprocessing raw data, including normalization and quality control.

Differential Modification Analysis: 
Code for identifying differentially unmodified positions (DUPs), differentially methylated positions (DMPs), and differentially hydroxymethylated positions (DHPs) across the genome.

Gene Ontology Enrichment: 
Scripts to perform gene ontology enrichment analysis for significant loci.

Validation: 
Code used for the technical validation of significant findings using bisulfite pyrosequencing data.

Single-Cell Analysis: 
Scripts for analyzing cell type-specific DNA modifications in serotonergic and non-serotonergic cells isolated from the DRN using laser capture microdissection (LCM).

Requirements

R version 3.3.2 or higher
Required R packages: wateRmelon, minfi, MLML2R, comb-p, missMethyl, among others.

Citation

Please cite the following publication when using this code:

Riemens, R.J.M., Pishva, E., Iatrou, A., et al. (2023). DNA modification profiling in the dorsal raphe nucleus highlights cell type-specific changes in TNXB in Alzheimer’s disease.bioRxiv. 2023:2023.08.28.555168.

For questions or issues, please contact e.pishva at e.pishva@maastrichtuniversity.nl.
