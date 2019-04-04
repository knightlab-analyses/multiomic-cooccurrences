## Michelle Badri
## 03/15/19
## Testing SE with microbe and metabolite data

setwd("~/depth-benchmarks")

library(phyloseq)
library(SpiecEasi)
library(igraph)


## Load converted tsv files, scrub any # or ' characters before import
metab <- read.table("table.rel.metabolites.1_a_scrub.txt", sep="\t",header=TRUE, row.names = 1)
microb <- read.table("table.rel.microbes.1_a_scrub.txt", sep="\t",header=TRUE, row.names = 1)

## Reorder the samples in one file to match the other
metab <- metab[,colnames(microb)]

## Load into phyloseq, object used to hold OTU tables
metabolite_phy <- phyloseq(OTU=otu_table(metab, taxa_are_rows=TRUE))
microbe_phy    <- phyloseq(OTU=otu_table(microb, taxa_are_rows=TRUE))

## Run SE on both objects
se.test <- spiec.easi(list(t(microbe_phy@.Data), t(metabolite_phy@.Data)), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05, ncores=2))

## Run SE on both objects
#se.test <- spiec.easi(list(t(microb), t(metab)), method='mb', nlambda=40,
#                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

## Plot using igraph 
dtype <- c(rep(1,ntaxa(microbe_phy)), rep(2,ntaxa(metabolite_phy)))
nodenames <- c(taxa_names(microbe_phy), taxa_names(metabolite_phy))
ig.se <- adj2igraph(getRefit(se.test))
lay = layout.fruchterman.reingold(ig.se)

#pdf("ig.se.test.depth.pdf",width = 10, height = 10)
plot(ig.se,layout=lay, vertex.color=dtype+1, 
     vertex.size=6, vertex.label=nodenames, vertex.label.cex=0.6,
    vertex.label.color="black",rescale=TRUE,vertex.label.family="Helvetica",
    ylim=c(-1,1),xlim=c(-1,1), asp = 0)
#dev.off()

#pdf("ig.se.test.nolabel.pdf",width = 10, height = 10)
plot(ig.se,layout=lay, vertex.color=dtype+1, 
     vertex.size=6, vertex.label=NA,
     vertex.label.color="black", rescale=TRUE,vertex.label.family="Helvetica",
     ylim=c(-1,1),xlim=c(-1,1), asp = 0)
#dev.off()

## Write edgelist
V(ig.se)$name <- nodenames
ig.el<- as_edgelist(ig.se,names=TRUE)
write.csv(ig.el,"edgelist_SEmultitest_depth.csv")


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
write.csv(ig.el.weight,"edgelist_weighted_SEmultitest_depth.csv")



