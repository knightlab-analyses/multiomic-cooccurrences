setwd("/Users/jmorton/Documents/dev/multiomic-cooccurences/ipynb")

library(propr)
library(magrittr)

met <- read.delim(
  "../data/soils/metabolites.tsv",
  skip = 1, row.names = 1) %>% t
mic <- read.delim(
  "../data/soils/microbes.tsv",
  skip = 1, row.names = 1) %>% t

# Fix row names
mic <- mic[rownames(met),]

# Replace with smallest non-zero value
met[met == 0] <- min(met[met != 0])
mic[mic == 0] <- min(mic[mic != 0])

# Do a CLR of met and mic separately
clr <- function(x) sweep(log(x), 1, rowMeans(log(x)), "-")
agg <- cbind(clr(met), clr(mic))

# PHI
pr <- propr:::lr2phi(agg)
colnames(pr) <- colnames(agg)
rownames(pr) <- colnames(agg)
write.csv(pr, "../results/soil_output/rerun/prop_matrix_soil_PHI.csv")

# RHO
pr <- propr:::lr2rho(agg)
colnames(pr) <- colnames(agg)
rownames(pr) <- colnames(agg)
write.csv(pr, "../results/soil_output/rerun/prop_matrix_soil_RHO.csv")


# PEARSON
pr <- stats::cor(agg)
colnames(pr) <- colnames(agg)
write.csv(pr, "../results/soil_output/rerun/prop_matrix_soil_pearson.csv")

# SPEARMAN
pr <- stats::cor(agg, method = "spearman")
colnames(pr) <- colnames(agg)
write.csv(pr, "../results/soil_output/rerun/prop_matrix_soil_spearman.csv")

