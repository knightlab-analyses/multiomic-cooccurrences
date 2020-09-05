# multiomic-cooccurences

This repository holds the analysis scripts and cystic fibrosis datasets to showcase the utility of applying neural networks to learn multi-omics cooccurences.

# Run scripts
The scripts within the `scripts` folder contains scripts used to run `rhapsody mmvec` for co-occurrrence analysis and `songbird multinomial` for differential abundance.
Case study specific visualizations and preprocessing can be found within the `ipynb` notebooks and in the `results`.

The command used to generate the cystic fibrosis co-occurences is found under the scripts folder. This was a result of multiple rounds of cross validation with parameters with learning rates (1e-5, 1e-6, 1e-7), input priors (0.1, 1, 10),  output priors (0.1, 1, 10) and a rank of 3.

When running these commands, it will take a substantial amount of compute time, so be sure to allocate multiple days to run.
It is also recommended to run Tensorboard in parallel in order to monitor the convergence of the algorithm.

# Notebooks
The notebooks for producing the figures in the manuscript.  They are given as follows
 - ipynb/Figure2-simulation.ipynb                  : Illustration of the CF biofilm simulation
 - ipynb/Figure2_biofilm.ipynb                     : CF biofilm simulation with addition sampling noise
 - ipynb/Figure2-benchmarks.ipynb                  : Benchmark figures based on CF biofilm simulation
 - ipynb/Figure3-biocrust.ipynb                    : Biocrust case study
 - ipynb/Figure4-cystic-fibrosis.ipynb             : Cystic Fibrosis case study
 - ipynb/Figure5-high-fat-diet.ipynb               : High fat diet case study
 - ipynb/Figure5-inflammatory-bowel-disease.ipynb  : IBD multiomics case study
 - ipynb/FigureS2-simulation.ipynb                 : Scaling simulation / benchmark

There are also auxiliary notebooks on how to extract parameters from a model checkpoint.
This is useful for real time diagnostics and debugging.

- ipynb/recover-CF-parameters.ipynb
- ipynb/biplot-coordinates.ipynb

# Results
This folder contains the checkpoints and diagnostics generated from the command provided in the scripts folder.

benchmark_output : simulation dataset used for the roc-curve analysis on mmvec, SPIEC-EASI, SparCC, proportionality,  pearson and spearman
soil_output      : results from the biocrust case study
cf_output        : results from the cystic fibrosis case study
hfd_output       : results from the high fat diet case study
ihmp_output      : results from the IBD case study

# Data files
All the data required to reproduce the analysis can be found under `data/`.

 - cf_sim : cystic fibrosis biofilm simulations generated from partial differential equations run in matlab using https://github.comoo/zhangzhongxun/WinCF_model_Code
 - CF     : cystic fibrosis data on oxygen gradients.  This includes 16S, metabolomics, taxonomies, sample metadata and LCMS annotations
 - soils  : 16S and metabolomics abundances for biocrust study
 - HFD    : high fat diet study on mice. 16S, metabolomics, taxonomies, sample metadata and LCMS annotations
 - ihmp   : IBD multiomicscase study. Taxonomic profiles from metagenomics, 4 metabolomics datasets, sample metadata and LCMS annotations
 
 # Rebuttal scripts
 All of the scripts used in our response are below
 https://github.com/knightlab-analyses/multiomic-cooccurrences/blob/rebuttal/ipynb/figure3-rerun.R
 https://github.com/knightlab-analyses/multiomic-cooccurrences/blob/rebuttal/ipynb/bayesian-optimization.ipynb
 https://github.com/knightlab-analyses/multiomic-cooccurrences/blob/rebuttal/ipynb/ground-truth-benchmark.ipynb
