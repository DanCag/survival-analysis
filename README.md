# Survival analysis
 In this github repository you will find data and scripts to reproduce the survival analysis of the [study](https://www.biorxiv.org/content/10.1101/2021.12.14.472046v1.full)
 
## Repository structure
 
`/data` folder:
- `/data/cBioPortal_survival-tables` contains the survival tables downloaded from [cBioPortal](https://www.cbioportal.org/) (TCGA, Firehose Legacy) for Colorectal Adenocarcinoma (COAD) and Lung Adenocarcinoma (LUAD) 
- `/data/CIBERSORTx` holds the inferred cell infiltration for TCGA-COAD and TCGA-LUAD patients. The cell infiltration is obtained by using the [CIBERSORTx](https://cibersortx.stanford.edu/) `Impute Cells Fractions` module. Two kinds of signature were used:
 - the LM22 signature that distinguishes 22 immun cell subtypes (provided by CIBERSORTx)
 - two signatures built after running CIBERSORTx module `Create Signature Matrix`. These signatures were built from the single-cell raw counts matrix obtained from sequencing and they were able to distinguish between 11 (_CD8-no-split_) and 12 (_CD8-split_) different immune cells types
- `/clinical` contains the clinical file downloaded from GDC for TCGA-COAD and TCGA-LUAD cohorts
- `/cytolytic-acitivty` contains the cytolytic activity of TCGA-COAD samples. Cytolytic activity was obtained from gene expression data (TPM) by computing the geometric mean of two genes: granzyme A (GZMA) and perforin (PRF1), as described in [Rooney et al., 2015](https://www.sciencedirect.com/science/article/pii/S0092867414016390?via%3Dihub)
- `/sample` contains the biospecimen file downloaded from [GDC](https://portal.gdc.cancer.gov/) for TCGA-COAD and TCGA-LUAD cohorts

<br>

`/output` folder: 
- processed survival tables with time and event information for TCGA-COAD and TCGA-LUAD (`surv-dfs_coad.rds`, `surv-os_luad.rds`). _OS_ stands for Overall Survival and the event is death, while _DFS_ stands for Disease-Free survival and the event is relapse
- processed survival tables to which we integrated molecular features such as cells infiltration information for each patient (_cells-fractions_) or both cells infiltration and cytolytic actvity info (_cells-fractions-cytolytic-activity_) for TCGA-COAD and TCGA-LUAD cohorts

<br>

`/R` folder:
1. `01_build_survival-dataframe.R` is a script that builds a survival dataframe starting from cBioPortal survival tables selecting for specific sample type and stage.
2. `02_integrate_cells-fractions.R` is a script that integrates cells infiltration to the survival dataframe built with 1.
3. `02_integrate_cells-fractions_cytolytic-activity.R` is a script that integrates cells infiltration and cytolytic activity to the survival dataframe built with 1.
4. `03_km.R` performs Kaplan-Meier analysis

<br>

## Software
 R version 4.2.0 (2022-04-22)  
 
 Packages:
 - `survival` (3.2-13)
 - `survminer` (0.4.9)
