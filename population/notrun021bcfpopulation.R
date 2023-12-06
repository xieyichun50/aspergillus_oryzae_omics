library(tidyr)
library(dplyr)
library(Hmisc)
library(stringr)
library(ggplot2)

setwd("202309aoryzae10/variant_bcf/")

samplelist<-read.delim("/mnt/content_176/yichun/fungi/aspergillus/202309aoryzae10/samples_n_reads_described.txt", 
                       header = F)
names(samplelist)<-c("group","sampleID")
vcf.merge<-as.data.frame(matrix(NA, nrow = 0, ncol = 4))
names(vcf.merge)<-c("CHROM","POS","REF","ALT")

for (i in 1:nrow(samplelist)){
  vcf.sub<-read.delim(paste0(samplelist$sampleID[i],".snp.vcf"),
                      header = F, comment.char = "#" )
  names(vcf.sub)<-c("CHROM","POS","ID",
                    "REF","ALT","QUAL","FILTER",
                    "INFO","FORMAT",
                    "detail")
  vcf.sub<-separate(vcf.sub, "INFO", c("DP", "I16"),
                    sep = ";", remove = TRUE, convert = TRUE)
  vcf.sub$DP<- gsub("DP=", "", vcf.sub$DP)
  vcf.sub$I16<-gsub("I16=", "", vcf.sub$I16)
  vcf.sub<-separate(vcf.sub, "I16", c("RF", "RR", "AF", "AR"),
                    sep = ",", remove = TRUE, convert = TRUE)
  #vcf.sub$DP<-as.numeric(vcf.sub$DP)
  vcf.sub$RF<-as.numeric(vcf.sub$RF)
  vcf.sub$RR<-as.numeric(vcf.sub$RR)
  vcf.sub$AF<-as.numeric(vcf.sub$AF)
  vcf.sub$AR<-as.numeric(vcf.sub$AR)
  vcf.sub$refr<-(vcf.sub$RF+vcf.sub$RR)/(vcf.sub$RF+vcf.sub$RR+vcf.sub$AF+vcf.sub$AR)
  #vcf.sub$pdepth<-vcf.sub$RF+vcf.sub$RR+vcf.sub$AF+vcf.sub$AR
  vcf.sub<-vcf.sub[(vcf.sub$RF+vcf.sub$RR+vcf.sub$AF+vcf.sub$AR)>=10 &
                     (vcf.sub$AF+vcf.sub$AR)>=5,
                   c("CHROM","POS","REF","ALT","refr")]
  names(vcf.sub)<-c("CHROM","POS","REF","ALT",samplelist$sampleID[i])
  vcf.merge<-merge(vcf.merge, vcf.sub, by = c("CHROM","POS","REF","ALT"), all = T)
}
rm(i, vcf.sub)
save.image("bcfmerge.RData")
write.table(vcf.merge, "aory.merge.reffq.txt",
            sep = "\t",
            quote = F, row.names = F)

##correlation
vcf.merge.matrix<-vcf.merge
row.names(vcf.merge.matrix)<-paste0(vcf.merge.matrix$CHROM,":",vcf.merge.matrix$POS,":",
                                    vcf.merge.matrix$REF,"|",vcf.merge.matrix$ALT)
vcf.merge.matrix<-vcf.merge.matrix[,5:ncol(vcf.merge.matrix)]

d<-cor(vcf.merge.matrix, method = "pearson", use = "na.or.complete")
View(as.data.frame(d))
h<-hclust(as.dist(1-d), method = "complete")

png(filename = paste0("aory.sample_pearson_tree.png"), 
    res = 300, units = "in", width = 4, height = 6)
plot(varclus(d, similarity="pearson", type = "similarity.matrix"), hang = -1)
dev.off()
rm(pca.eigenvalues,pca.matrix,pca.scores,h,p,d)
