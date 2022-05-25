# Survival analysis
 In this github repository you will find data and scripts to reproduce the survival analysis of the [study](https://www.biorxiv.org/content/10.1101/2021.12.14.472046v1.full)
 
## Repository structure
 
`/data` folder:
- `/data/cBioPortal_survival-tables` contains the survival tables downloaded from [cBioPortal](https://www.cbioportal.org/) (TCGA, Firehose Legacy) for Colorectal Adenocarcinoma (COAD) and Lung Adenocarcinoma (LUAD) 
- `/data/CIBERSORTx` holds the inferred cell infiltration for TCGA-COAD and TCGA-LUAD patients. The cell infiltration is obtained by using the [CIBERSORTx](https://cibersortx.stanford.edu/) `Impute Cells Fractions` module. Two kinds of signature were used:
 - the _LM22_ signature that distinguishes 22 immune cell subtypes (provided by CIBERSORTx)
 - two signatures built after running CIBERSORTx module `Create Signature Matrix`. These signatures were built from the raw counts matrix obtained from single-cell sequencing and they were able to distinguish between 11 (_CD8-no-split_) and 12 (_CD8-split_) different immune cells types
- `/clinical` contains the clinical file downloaded from [GDC](https://portal.gdc.cancer.gov/) for TCGA-COAD and TCGA-LUAD cohorts
- `/cytolytic-acitivty` contains the cytolytic activity of TCGA-COAD samples. Cytolytic activity was obtained from gene expression data (TPM) by computing the geometric mean of two genes: granzyme A (GZMA) and perforin (PRF1), as described in [Rooney et al., 2015](https://www.sciencedirect.com/science/article/pii/S0092867414016390?via%3Dihub)
- `/sample` contains the biospecimen file downloaded from [GDC](https://portal.gdc.cancer.gov/) for TCGA-COAD and TCGA-LUAD cohorts

<br>

`/R` folder contains scripts:
- `01_build_survival-dataframe.R`: builds a survival dataframe starting from cBioPortal survival tables selecting for specific sample type and stage.
- `02_integrate_cells-fractions.R`: integrates cells infiltration to the survival dataframe built with `01_` employing either the _CD8-split_ or _CD8-no-split_ signature.
- `02_integrate_cells-fractions_cytolytic-activity.R` integrates cells infiltration (employing the _LM22_ signature) and cytolytic activity to the survival dataframe built with `01_`.
- `03_km.R` performs Kaplan-Meier analysis

<br>

`/output` is the folder where you will store the ouptut files generated after running scripts `01_` and `02_`.<br>
In particular you will generate:
- processed survival tables with time and event information for TCGA-COAD and TCGA-LUAD (`surv-dfs_coad.rds`, `surv-os_luad.rds`) where _OS_ stands for _Overall Survival_ (death is the event), while _DFS_ stands for _Disease-Free survival_ (relapse is the event). These are the output of `01_build_survival-dataframe.R`
- processed survival tables to which molecular features such as cells infiltration information for each patient (_cells-fractions_) or both cells infiltration and cytolytic actvity info (_cells-fractions-cytolytic-activity_) are integrated for TCGA-COAD and TCGA-LUAD cohorts. These are the output of `02_integrate_cells-fractions.R` and `02_integrate_cells-fractions_cytolytic-activity.R` respectively.
- `.Rdata` files with parameters saved from script `01_` (`build_parameters.RData`) and from script `02_` (`integrate_parameters.RData`)

<br>

`/example_output` folder contains the output files (survival table + cells fractions or survival table + cells fraction-cytolytic activity info) you generate after running scripts `01_` and `02_`. You can use these files to check if everything is going smooth while running the scripts.



<br>


## Instructions
The numbers in the scripts' names indicate in which order they need to be run. <br>

For the second step there are two alternative scripts: `02_integrate_cells-fractions.R` integrates cells infiltration on the processed survival table, while `02_integrate_cells-fractions_cytolytic-activity.R` integrates both cells infiltration and cytolytic activity to the processed survival table. Depending on which analysis you want to reproduce you will need to run one or the other `02_` script with specific parameters (see _Reproduce plots_ section).

### Reproduce plots
- Figure 6G and Figure Supp. 6I<br>

Parameters of `01_build_survival-dataframe.R` script are 
```
tumor_type <- "COAD"
survival_analysis <- "DFS"
stages_info <- "I_IV"
```
In this case, you use `02_integrate_cells-fractions.R` script with following parameters
```
signature <- "CD8-split"
feature <- "cells-fractions"
```

Parameters of `03_km.R` script is
```
strat <- "CD8_Tem_GZMK_high"
```
Please note that to plot Figure Supp. 6I you need to uncomment and run line 130 in `03_km.R` script.


- Figure Supp. 2A<br>

Paramaters of `01_build_survival-dataframe.R` are
```
tumor_type <- "COAD"
survival_analysis <- "DFS"
stages_info <- "I_IV"
```

In this case you use `02_integrate_cells-fractions_cytolytic-activity.R` script with following parameters
```
signature <- "LM22"
feature <- "cells-fractions_cytolytic-activity"
```

Parameters of `03_km.R` script is
```
strat <- "Cyt"
```

- Figure Supp. 4A<br>

Paramaters of `01_build_survival-dataframe.R` are
```
tumor_type <- "COAD"
survival_analysis <- "DFS"
stages_info <- "I_IV"
```

In this case you use `02_integrate_cells-fractions_cytolytic-activity.R` script with following parameters
```
signature <- "LM22"
feature <- "cells-fractions_cytolytic-activity"
```

Parameters of `03_km.R` script is
```
strat <- "Neutrophils_Cyt"
```

- Figure Supp. 6J<br>

Parameters of `01_build_survival-dataframe.R` script are 
```
tumor_type <- "LUAD"
survival_analysis <- "OS"
stages_info <- "I_IV"
```
In this case you use `02_integrate_cells-fractions.R` script with following parameters
```
signature <- "CD8-split"
feature <- "cells-fractions"
```

Parameters of `03_km.R` script is
```
strat <- "CD8_Tem_GZMK_high"
```

<br>

## Software
 R version 4.2.0 (2022-04-22)  
 
 Packages:
 - `survival` (3.2-13)
 - `survminer` (0.4.9)

<br>

## Operating system
Ubuntu Ubuntu 20.04.4 LTS
