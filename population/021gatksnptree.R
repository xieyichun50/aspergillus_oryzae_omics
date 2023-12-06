library(gdsfmt)
library(SNPRelate) 
library(ggplot2)
library(ggtree)
library(ape)

setwd("202309aoryzae10/variant/")
##read in vcf (including koji ABCD)
snpgdsVCF2GDS("aory.merge.vcf",
              "aory.merge.gds",
              method = "biallelic.only")
input.gds<- snpgdsOpen("aory.merge.gds")
input.ibs<- snpgdsHCluster(snpgdsIBS(input.gds,num.thread=2, autosome.only=FALSE))
input.cuttree<-snpgdsCutTree(input.ibs)
treefile<- input.cuttree$dendrogram

png(filename = paste0("gatksnptree.png"), 
    res = 300, units = "in", width = 6, height = 4)
plot(input.cuttree$dendrogram,horiz=T)
dev.off()

##read in vcf
snpgdsVCF2GDS("aory.2023.vcf",
              "aory.2023.gds",
              method = "biallelic.only")
input.gds<- snpgdsOpen("aory.2023.gds")
input.ibs<- snpgdsHCluster(snpgdsIBS(input.gds,num.thread=2, autosome.only=FALSE))
input.cuttree<-snpgdsCutTree(input.ibs)
treefile<- input.cuttree$dendrogram

png(filename = paste0("gatksnptree.2023.png"), 
    res = 300, units = "in", width = 6, height = 4)
plot(input.cuttree$dendrogram,horiz=T)
dev.off()
