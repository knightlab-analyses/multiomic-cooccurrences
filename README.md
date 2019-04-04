# multiomic-cooccurences

This repository holds the analysis scripts and cystic fibrosis datasets to showcase the utility of applying neural networks to learn multi-omics cooccurences.

# Run scripts



The command used to generate the cystic fibrosis occurences is found under the scripts folder. This was a result of multiple rounds of cross validation with parameters with learning rates (1e-5, 1e-6, 1e-7), input priors (0.1, 1, 10),  output priors (0.1, 1, 10) and a rank of 3.

When running these commands, it will take a substantial amount of compute time, so be sure to allocate multiple days to run.
It is also recommended to run Tensorboard in parallel in order to monitor the convergence of the algorithm.

# Notebooks
The notebooks for producing the figures in the manuscript.  They are given as follows
 - ipynb/Figure2/biofilm.ipynb: The creation of the biofilm simulation
 - ipynb/Figure2/roc-curve.ipynb: The roc curve analysis
 - ipynb/Figure2/scale-benchmarks.ipynb: The scale-invariance benchmark analysis
 - ipynb/recover-CF-parameters.ipynb : Extract the parameters from the model.
 - ipynb/biplot-coordinates.ipynb : Calculate biplot coordinates to visualize co-occurence probabilities in Emperor
 - ipynb/Figure3.ipynb : The cystic fibrosis study in Figure 3
 - ipynb/Figure4.ipynb : The biocrust wetting experiment in Figure 4
 - ipynb/Figure5.ipynb : The high fat diet study in Figure 5

# Results
This folder contains the checkpoints and diagnostics generated from the command provided in the scripts folder.

 - depth_benchmark : simulation dataset used for the roc-curve analysis on mmvec, spiec-easi, pearson and spearman
 - effect_size_benchmark : simulation dataset used for the scale invariance analysis on mmvec, spiec-easi, pearson and spearman
 - cf_output : cystic fibrosis model results with mmvec and pearson
 - soil_output : biocrust model results with mmvec, spearman and spiec-easi
 - hfd_output : high fat diet model results with mmvec and pearson

# Data files
All the data required to reproduce the analysis can be found under `data/`.

 - cf_sim : cystic fibrosis biofilm simulations generated from partial differential equations run in matlab using https://github.com/zhangzhongxun/WinCF_model_Code
 - CF : cystic fibrosis data on oxygen gradients.  This includes 16S, metabolomics, taxonomies, sample metadata and LCMS annotations
 - soils : 16S and metabolomics abundances for biocrust study
 - HFD : high fat diet study on mice. 16S, metabolomics, taxonomies, sample metadata and LCMS annotations
