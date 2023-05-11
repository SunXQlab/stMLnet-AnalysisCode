# Modeling spatial intercellular communication and multilayer signaling regulations using stMLnet

## Introduction
Multicellular organisms require intercellular and intracellular signaling to coordinately regulate different cell functions. The technological advance of spatial transcriptomics (ST) lets us leverage spatial information to better elucidate cell signaling and functioning. 

Here, we present stMLnet, a method that infers spatial intercellular communication and multilayer signaling regulations from ST data by quantifying distance-weighted ligand-receptor (LR) signaling activity based on diffusion and mass action models and mapping it to intracellular targets. 

We demonstrated the applicability of stMLnet on multiple datasets, such as breast cancer and COVID-19 microenvironment. We enchmarked stMLnet's performance using multiple cell line perturbation datasets, synthetic data, and LR-target correlations stratified by cellular distance. Furthermore, we applied stMLnet to a scRNA-seq dataset of gliomas for trying to expand the scope of stMLnet.

The R package of stMLnet is publicly available from <a href="https://github.com/SunXQlab/stMLnet" target="_blank">stMLnet repository</a>. A web-based application of stMLnet is also developed and available at www.stmlnet.top/net. 
 
## Workflow

1.  **apply_in_stBC** contains the code to reproduce plots and benchmarking of the breast cancer dataset <br>
    - s1_runMLnet.R: construct a multilayer signaling network in breast cancer environment.
    - s2_calculate_LRTG_activity.R: infer LR signling activate based on the expression and distance of ligands and receptor.
    - s3_calculate_LRTG_importance.R: train a random forest model to predicte LR-target gene regulation.
    - s4_compare_method.R: compare the performance of different distance algorithm on breast cancer datasets, corresponding to Fig3D.
    - s5_compare_software.R: compare the performance of similar software (MISTy, NicheNet, CytoTalk) on breast cancer datasets, corresponding to Fig2F.
    - s6_visualize_CCI.R: various visualizations of the commucation in breast cancer environment, corresponding to Fig2B-E, FigS3.
    - s7_get_cor.R: calculate the LR-target correlations,responding to FigS4.
2. **apply_in_simu** contains the code to reproduce the simulation study of stMLnet, corresponding to Fig3C. <br>
3. **apply_in_scST** contains the code to reproduce the plot and detailed analysis of the three single-cell resolution ST datasets.<br>
    - **giotto_seqfish_dataset** contains the code to reproduce the plot and detailed analysis on seqfish+ dataset.<br>
        + s4_visial_CCI.R: various visualizations of the commucation in seqfish+ dataset, corresponding to Fig4B-E.
    - **giotto_slideseq2_dataset** contains the code to reproduce the plot and detailed analysis on slide-seq2 dataset.<br>
        + s4_visial_CCI.R: various visualizations of the commucation in slide-seq2 dataset, corresponding to Fig5B-C.
   - **giotto_merfish_dataset** contains the code to reproduce the plot and detailed analysis on merfish dataset.<br>
        + s0_compare_merfish_results.R: comparison of cell communication in different layers of MERFISH data, corresponding to FigS6.
        + giotto_merfish_dataset_layer9/s4_visial_CCI.R: various visualizations of the commucation in merfish dataset, corresponding to Fig5F-H. 
4. **apply_in_stGBM** contains the code to reproduce plots and detailed analysis for appling stMLnet on the gliomas dataset <br>
   - s4_visial_CCI.R: various visualizations of the commucation in gliomas environment, corresponding to Fig6B-C, FigS7.
   - s5_compare_software.R: compare the performance of similar software (NicheNet, CytoTalk) on gliomas datasets, corresponding to Fig6E.
   - s6_gesaAnalysis.R: GSEA analysis of macrophages and malignnat, corresponding to FigS8.
5. **apply_in_COVID19** contains the code to reproduce plots and detailed analysis of the COVID-19 ST dataset <br>
   - s4_visial_CCI.R: various visualizations of the commucation in breast cancer environment, corresponding to Fig7B-C, Fig7F, FigS10.
   - s5_check_feedback_loop.R: visualizations of the positive feedback circuits between AECs, macrophages and monocytes, corresponding to Fig7D, FigS11.
6. **code** contains all functions of stMLnet to analysis cell-cell interactions <br>

## Databases

Databases used in this paper are stored in `./prior_knowledge/output` folder, including LigRecDB (the Ligand-Receptor database), RecTFDB (Receptor-TF database) and TFTGDB (TF-Target Gens database). The R code used for collection and integration of prior databases is available at `./prior_knowledge/`.

  - s1_get_prior_database.R: collection and integration of databases, corresponding to FigS1.
  - s2_compare_database.R: comparison of different databases (stMLnet, Omnipath, NicheNet), corresponding to FigS2A.
  - s3_tunepara_database.R: train the hyperparameters of random walk with restart on different databases, corresponding to FigS2B.
  - s4_doRWR.R: perform the random walk with restart to obtain the corrected databases.
  - s5_quan.cutoff.R: Using the cell lines databases to optimizate the corrected databases, corresponding to FigS2C.
    

