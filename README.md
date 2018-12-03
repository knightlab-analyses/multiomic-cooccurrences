# multiomic-cooccurences

This repository holds the analysis scripts and cystic fibrosis datasets to showcase the utility of applying neural networks to learn multi-omics cooccurences.

# Run scripts

The command used to generate the cystic fibrosis occurences is found under the scripts folder. This was a result of multiple rounds of cross validation with parameters with learning rates (1e-5, 1e-6, 1e-7), input priors (0.1, 1, 10),  output priors (0.1, 1, 10) and a rank of 3.

When running these commands, it will take a substantial amount of compute time, so be sure to allocate multiple days to run.
It is also recommended to run Tensorboard in parallel in order to monitor the convergence of the algorithm.

# Notebooks
The notebooks are mainly for analyzing the fitted model.  They are given in the following order
- recover-parameters.ipynb : Extract the parameters from the model.  Note that this can be performed in the middle of performing a fit
- biplot-coordinates.ipynb : Calculate biplot coordinates to visualize co-occurence probabilities in Emperor
- cf-figures.ipynb : Recreate the figures shown in the original paper


# Data files
- lcms_annotations.txt : GNPS annotations for MS2 spectra
- lcms_nt.biom : MS1 bucket table of metabolite abundances
- otus_nt.biom : Deblurred biom for microbe abundances
- taxonomy.tsv : taxonomic annotations for microbes
- validated-molecules.csv : validated annotations for spectra