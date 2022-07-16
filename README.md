# Modeling spatial intercellular communication and multilayer signaling regulations using stMLnet

## Introduction 
Multicellular organisms require intercellular and intracellular signaling to coordinately regulate different cell functions. The technological advance of spatial transcriptomics (ST) lets us leverage spatial information to better elucidate cell signaling and functioning. 

Here, we present stMLnet, a method that infers spatial intercellular communication and multilayer signaling regulations from ST data by quantifying distance-weighted ligand-receptor (LR) signaling activity based on diffusion and mass action models and mapping it to intracellular targets. 

We demonstrated the applicability of stMLnet on multiple datasets, such as breast cancer and COVID-19 microenvironment. We enchmarked stMLnet's performance using multiple cell line perturbation datasets, synthetic data, and LR-target correlations stratified by cellular distance. Furthermore, we applied stMLnet to a scRNA-seq dataset of gliomas for trying to expand the scope of stMLnet.

The R package of stMLnet is publicly available from <a href="https://github.com/SunXQlab/stMLnet" target="_blank">stMLnet repository</a>. A web-based application of stMLnet is also developed and available at www.stmlnet.top/net. 
 
## Workflow

1. **apply_in_stBC** contains the code to reproduce plots and benchmarking of the breast cancer dataset <br>
2. **apply_in_COVID19** contains the code to reproduce plots and detailed analysis of the COVID-19 ST dataset <br>
3. **apply_in_simu** contains the code to reproduce the simulation study of stMLnet <br>
4. **apply_in_scGBM** contains the code to reproduce plots and detailed analysis for appling stMLnet on the scRNA-seq dataset <br>
5. **code** contains all functions of stMLnet to analysis cell-cell interactions <br>

## Databases

Databases used in this paper are stored in `./prior_knowledge/output` folder, including LigRecDB (the Ligand-Receptor database), RecTFDB (Receptor-TF database) and TFTGDB (TF-Target Gens database). The R code used for collection and integration of prior databases is available at `./prior_knowledge/`.
