# FMT-IBD-project

This repository contains the reproducible analysis for the study:

**Dynamic Alterations in Stool Microbiota Composition Following Fecal Transplantation in Patients with Inflammatory Bowel Disease**

![](https://github.com/JeniaOle13/FMT_IBD_project/blob/main/manuscript/Figure_1.png)
Study design (Created with BioRender.com). (A) Schematic representation of the patient treatment scheme; colours denote disease groups. (B) Longitudinal stool sampling scheme; the magenta marker indicates the time point of the FMT procedure. (C) Stool sampling time points for each patient, with patient IDs shown on the left; colours denote disease groups. FMT, fecal microbiota transplantation.

Quarto report available [here](https://jeniaole13.github.io/FMT_IBD_project/)

The project investigates the effects of fecal microbiota transplantation (FMT) on the gut microbiota of patients with inflammatory bowel disease (IBD), focusing on differences between Crohn's disease (CD) and ulcerative colitis (UC) within one month after FMT.

## What is included

- A **Quarto document** (`main.qmd`) that reproduces all analyses, from taxonomic and functional profiles processing to figures, tables, and statistical tests.
- Code for taxonomic and functional profiling, alpha and beta diversity, DEICODE-based trajectory analysis, Songbird differential ranking, and RECAST donor-derived species detection.
- Scripts for statistical modeling (mixed-effects and generalized linear models) and post-hoc power analysis.
- Intermediate and supplementary data tables used in the manuscript.

## Key findings

- UC patients showed a coherent shift of the microbiota towards the donor and towards healthy controls, with significant temporal changes in taxonomic and functional profiles.
- CD patients showed less pronounced taxonomic restructuring despite comparable clinical improvement, with only a weak functional signal.
- Power analysis indicated that the RPCA-based interaction was adequately powered, whereas the time × disease interactions for distance from donor and colonisation were underpowered, limiting conclusions for those endpoints.

## Data and tools

- Raw metagenomic reads are available in the NCBI repository under project ID [PRJNA763503](https://www.ncbi.nlm.nih.gov/bioproject/763503).
- The analysis uses the [ASSNAKE](https://github.com/ASSNAKE) pipeline for metagenomic processing and [RECAST](https://github.com/ctlab/recast) for donor-derived species identification.
- All analyses were performed in R and reported using Quarto.

## Contact

For questions, please contact [Evgenii I. Olekhnovich](mailto:jeniaole01@gmail.com).
