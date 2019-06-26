setwd("~/Desktop/SEMultiomicsTest/reviseddepth_benchmarks/CF_sim")

## RELATIVE
library(phyloseq)
library(igraph)
library(devtools)
library(SpiecEasi)
library(propr)

metab  <- read.table("rel.metab.tsv",sep="\t",header=T,row.names = 1,quote="")
microb <- read.table("rel.microbes.tsv",sep="\t",header=T,row.names = 1)

## Reorder the samples in one file to match the other
microb <- microb[,colnames(metab)]


## Load into phyloseq, object used to hold OTU tables
metabolite_phy <- phyloseq(OTU=otu_table(metab, taxa_are_rows=TRUE))
microbe_phy    <- phyloseq(OTU=otu_table(microb, taxa_are_rows=TRUE))


## Run SE on both objects
se.test <- spiec.easi(list(t(microbe_phy@.Data), t(metabolite_phy@.Data)), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=2))


## Plot using igraph 
dtype <- c(rep(1,ntaxa(microbe_phy)), rep(2,ntaxa(metabolite_phy)))
nodenames <- c(taxa_names(microbe_phy), taxa_names(metabolite_phy))
ig.se <- adj2igraph(getRefit(se.test))

## Get edgelist weights from beta (optimal covariance) matrix 
sebeta <- as.matrix(symBeta(getOptBeta(se.test), mode='maxabs'))
rownames(sebeta) <- colnames(sebeta) <- nodenames
el <- get.edgelist(ig.se)
sizes <- list(rep(1, length(el[,1])))
for (i in 1:length(el[,1])){
  first <- el[,1][i]
  second <- el[,2][i] 
  sizes[i]<-sebeta[first,second]
}

E(ig.se)$weight <- unlist(sizes)


## Write edgelist with weights from beta matrix: model coefficients from neighborhood selection
V(ig.se)$name <- nodenames
ig.el<- as_edgelist(ig.se,names=TRUE)
ig.el.weight <- cbind(ig.el , round(E(ig.se)$weight,5 ))
write.csv(ig.el.weight,"edgelist_rel_0622_SEmultitest.csv")

sebeta <- as.matrix(symBeta(getOptBeta(se.test), mode='maxabs'))
rownames(sebeta) <- colnames(sebeta) <- nodenames
write.csv(sebeta,"weights_rel_0622_SPIECEASI.csv")
dim(metab)
dim(microb)

comb<- rbind(metab,microb)
sparc_test<- sparcc(t(comb),iter=10, th=0)
cor_result<- sparc_test$Cor
rownames(cor_result) <- colnames(cor_result) <- rownames((comb))

write.csv(cor_result, "cor_matrix_rel_fixed_0622_SparCC.csv")

comb<- rbind(metab,microb)
rho           <- propr(t(comb), metric = "rho")
rhoprop  <- rho@matrix

write.csv(rhoprop, "prop_matrix_rel_fixed_0622_RHO.csv")

phi           <- propr(t(comb), metric = "phi")
phiprop  <- phi@matrix

write.csv(phiprop, "prop_matrix_rel_fixed_0622_PHI.csv")





## ABSOLUTE

library(phyloseq)
library(igraph)
library(devtools)
library(SpiecEasi)
library(propr)

metab  <- read.table("abs.metab.tsv",sep="\t",header=T,row.names = 1,quote="")
microb <- read.table("abs.microbes.tsv",sep="\t",header=T,row.names = 1)

## Reorder the samples in one file to match the other
microb <- microb[,colnames(metab)]


## Load into phyloseq, object used to hold OTU tables
metabolite_phy <- phyloseq(OTU=otu_table(metab, taxa_are_rows=TRUE))
microbe_phy    <- phyloseq(OTU=otu_table(microb, taxa_are_rows=TRUE))


## Run SE on both objects
se.test <- spiec.easi(list(t(microbe_phy@.Data), t(metabolite_phy@.Data)), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=2))


## Plot using igraph 
dtype <- c(rep(1,ntaxa(microbe_phy)), rep(2,ntaxa(metabolite_phy)))
nodenames <- c(taxa_names(microbe_phy), taxa_names(metabolite_phy))
ig.se <- adj2igraph(getRefit(se.test))

## Get edgelist weights from beta (optimal covariance) matrix 
sebeta <- as.matrix(symBeta(getOptBeta(se.test), mode='maxabs'))
rownames(sebeta) <- colnames(sebeta) <- nodenames
el <- get.edgelist(ig.se)
sizes <- list(rep(1, length(el[,1])))
for (i in 1:length(el[,1])){
  first <- el[,1][i]
  second <- el[,2][i] 
  sizes[i]<-sebeta[first,second]
}

E(ig.se)$weight <- unlist(sizes)


## Write edgelist with weights from beta matrix: model coefficients from neighborhood selection
V(ig.se)$name <- nodenames
ig.el<- as_edgelist(ig.se,names=TRUE)
ig.el.weight <- cbind(ig.el , round(E(ig.se)$weight,5 ))
write.csv(ig.el.weight,"edgelist_abs_0622_SEmultitest.csv")

sebeta <- as.matrix(symBeta(getOptBeta(se.test), mode='maxabs'))
rownames(sebeta) <- colnames(sebeta) <- nodenames
write.csv(sebeta,"weights_abs_0622_SPIECEASI.csv")
dim(metab)
dim(microb)

comb<- rbind(metab,microb)
sparc_test<- sparcc(t(comb),iter=10, th=0)
cor_result<- sparc_test$Cor
rownames(cor_result) <- colnames(cor_result) <- rownames((comb))

write.csv(cor_result, "cor_matrix_abs_fixed_SparCC.csv")

comb<- rbind(metab,microb)
rho           <- propr(t(comb), metric = "rho")
rhoprop  <- rho@matrix

write.csv(rhoprop, "prop_matrix_abs_fixed_2_RHO.csv")

phi           <- propr(t(comb), metric = "phi")
phiprop  <- phi@matrix

write.csv(phiprop, "prop_matrix_abs_fixed_2_PHI.csv")


