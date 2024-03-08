# Modeling spatial intercellular communication and multilayer signaling regulations using stMLnet

## Introduction
Multicellular organisms require intercellular and intracellular signaling to coordinately regulate different cell functions. The technological advance of spatial transcriptomics (ST) lets us leverage spatial information to better elucidate cell signaling and functioning. 

We propose an ST-based multilayer network method, stMLnet, for inferring spatial intercellular communication and multilayer signaling regulations by quantifying distance-weighted ligand–receptor signaling activity based on diffusion and mass action models and mapping it to intracellular targets. To benchmark stMLnet, we designed a benchmark framework to evaluate and compare its performance with other representative CCC inference methods for inferring intercellular and intracellular communication. We compared stMLnet with six state-of-the-art CCC inference methods, including CellChatV2, COMMOT, CytoTalk, MISTy, NicheNet, and Scriabin. We evaluated these seven methods using three ST datasets and used 15 sets of cell line perturbed expression data as the ground truth to assess the prediction results. 

Furthermore, we demonstrate the applicability of stMLnet on six ST datasets acquired with four different technologies (e.g., seqFISH+, Slide-seq v2, MERFISH, Stereo-seq, and Visium), showing its effectiveness and reliability on ST data with varying spatial resolutions and gene coverages. Finally, stMLnet identifies positive feedback circuits between alveolar epithelial cells, macrophages, and monocytes via multilayer signaling pathways within a COVID-19 microenvironment. Our proposed method provides an effective tool for predicting multilayer signaling regulations between interacting cells, which can advance the mechanistic and functional understanding of spatial CCCs.

The R package of stMLnet is publicly available from <a href="https://github.com/SunXQlab/stMLnet" target="_blank">stMLnet repository</a>. 
 
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
        + s4_visial_CCI.R: various visualizations of the commucation in seqfish+ dataset, corresponding to Fig5A-C.
   - **giotto_merfish_dataset** contains the code to reproduce the plot and detailed analysis on merfish dataset.<br>
        + s0_compare_merfish_results.R: comparison of cell communication in different layers of MERFISH data, corresponding to FigS6.
        + giotto_merfish_dataset_layer6/s2_visial_CCI.R & s3_visial_spatialCCC.R: various visualizations of the communication in merfish dataset, corresponding to Fig5F and Fig5G-J, respectively.
   - **giotto_slideseq v2_dataset** contains the code to reproduce the plot and detailed analysis on slide-seq v2 dataset.<br>
        + s4_visial_CCI.R: various visualizations of the commucation in slide-seq v2 dataset, corresponding to Fig6B-F.
   - **giotto_Stereoseq_dataset** contains the code to reproduce the plot and detailed analysis on Stereo-seq dataset.<br>
4. **apply_in_stGBM** contains the code to reproduce plots and detailed analysis for appling stMLnet on the gliomas dataset <br>
   - s4_visial_CCI.R: various visualizations of the commucation in gliomas environment, corresponding to Fig6B-C, FigS7.
   - s5_compare_software.R: compare the performance of similar software (NicheNet, CytoTalk, MISTy) on gliomas datasets, corresponding to Fig6E.
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
    
## Reproducible run

We provide a reproducible run capsule in <a href="https://codeocean.com/capsule/9121262/tree" target="_blank">CodeOcean</a>, showing the application of stMLnet in SlideSeqV2 dataset. This capsule includes enviorment installation, execuate code and necessary data, taking about 1.35 hours to set up enviorment and to reproduce related results. (Note: we made minor changes of input and output path in code to adapt in docker images)

You can reproduce the paper results following the codes, please make sure you have already installed related dependencies:

       # Check if the following dependencies are installed.
       pkgs <- c('Seurat','SeuratWrappers','Giotto','reshape2','stringr','dplyr', # for data preprocessing
                        'caret','doParallel','snow','foreach', # for quantitative model
                         'ggplot2','ggsci','clusterProfiler','org.Hs.eg.db', 'plotrix','ggalluvial','ggraph','igraph' # for visualization
                         )
       for (pkg in pkgs) {
         if (!requireNamespace(pkg)) { cat(paste0('please install and library the package: ',pkg,'\n')) }
       }
       
       # Installing related dependencies.
       pkgs <- c( 'caret','doParallel','snow','foreach','ggplot2','ggsci','clusterProfiler','org.Hs.eg.db','plotrix','ggalluvial','ggraph','igraph')
       for (pkg in pkgs) {install.packages(pkg, repos = 'https://cloud.r-project.org')}
       
       devtools::install_version("spatstat.core", version = "2.4-4", repos="https://cloud.r-project.org/")
       devtools::install_version("Seurat", version = "4.2.0", repos="https://cloud.r-project.org/")
       remotes::install_github("satijalab/seurat-wrappers")
       remotes::install_github("drieslab/Giotto",  ref="v1.1.0")

If you have problems installing the environment manually, you can also choose to install the dependent environment via dockfile:

       # Bash
       # built a docker image
       # ensure that dockerfile and postInstall are in the same path
       docker bulid -f Dockerfile -t stMLnetEnv:0.1 .
       # Run docker image
       docker run -it stMLnetEnv:0.1 /bin/bash
