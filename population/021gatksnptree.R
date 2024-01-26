library(dplyr)
library(tidyr)
library(gdsfmt)
library(SNPRelate) 
library(pheatmap)
library(ggplot2)
library(ggtree)
library(ape)

setwd("202309aoryzae10/variant/")
####all samples + RIB40 (including kojiABCD)
samplegroup=".all_out"
sample.order<-c("TSS190621","TS3230418","AoL",
                "koji20230527","koji20230603",
                "koji20231118","koji20231209",
                "koji2018A","koji2018B","koji2018C","koji2018D",
                "AoH","RDS170320","RD3230418",
                "RIB40")
snpgdsVCF2GDS(paste0("aory",samplegroup,".vcf"),
              paste0("aory",samplegroup,".gds"),
              method = "biallelic.only")
input.gds<- snpgdsOpen(paste0("aory",samplegroup,".gds"))
input.diss<-snpgdsDiss(input.gds,
                       remove.monosnp=F,
                       autosome.only=F)

input.ibs<- snpgdsHCluster(snpgdsIBS(input.gds,num.thread=2, autosome.only=FALSE))
input.cuttree<-snpgdsCutTree(input.ibs)
treefile<- input.cuttree$dendrogram

par(mar = c(5, 5, 5, 5), xpd =T)
png(filename = paste0("aory",samplegroup,".gatksnptree.png"), 
    res = 300, units = "in", width = 6, height = 6)
plot(input.cuttree$dendrogram, horiz=T, cex = 0.8, xlim = c(0.005,-0.005))
dev.off()

####all samples (including kojiABCD)
samplegroup=".all"
sample.order<-c("TSS190621","TS3230418","AoL",
                "koji20230527","koji20230603",
                "koji20231118","koji20231209",
                "koji2018A","koji2018B","koji2018C","koji2018D",
                "AoH","RDS170320","RD3230418")
####only 2023 samples
samplegroup=".new_out"
sample.order<-c("TSS190621","TS3230418","AoL",
                "koji20230527","koji20230603",
                "koji20231118","koji20231209",
                "RDS170320","AoH","RD3230418",
                "RIB40")
snpgdsVCF2GDS(paste0("aory",samplegroup,".vcf"),
              paste0("aory",samplegroup,".gds"),
              method = "biallelic.only")
input.gds<- snpgdsOpen(paste0("aory",samplegroup,".gds"))
input.diss<-snpgdsDiss(input.gds,
                       remove.monosnp=F,
                       autosome.only=F)

input.ibs<- snpgdsHCluster(snpgdsIBS(input.gds,num.thread=2, autosome.only=FALSE))
input.cuttree<-snpgdsCutTree(input.ibs)
treefile<- input.cuttree$dendrogram

par(mar = c(5, 5, 5, 5), xpd =T)
png(filename = paste0("aory",samplegroup,".gatksnptree.png"), 
    res = 300, units = "in", width = 6, height = 6)
plot(input.cuttree$dendrogram, horiz=T, cex = 0.8, xlim = c(0.005,-0.005))
dev.off()

####only 2023 samples
samplegroup=".new"
sample.order<-c("AoH","RD3230418","RDS170320",
                "koji20230527","koji20230603",
                "koji20231118","koji20231209",
                "TSS190621","TS3230418","AoL")

##read in vcf (including koji ABCD)
snpgdsVCF2GDS(paste0("aory",samplegroup,".vcf"),
              paste0("aory",samplegroup,".gds"),
              method = "biallelic.only")
input.gds<- snpgdsOpen(paste0("aory",samplegroup,".gds"))
input.diss<-snpgdsDiss(input.gds,
                       remove.monosnp=F,
                       autosome.only=F)

input.ibs<- snpgdsHCluster(snpgdsIBS(input.gds,num.thread=2, autosome.only=FALSE))
input.cuttree<-snpgdsCutTree(input.ibs)
treefile<- input.cuttree$dendrogram

par(mar = c(5, 5, 5, 5), xpd =T)
png(filename = paste0("aory",samplegroup,".gatksnptree.png"), 
    res = 300, units = "in", width = 6, height = 4)
plot(input.cuttree$dendrogram, horiz=T, cex = 0.8, xlim = c(0.4,-0.05))
dev.off()

##similarity matrix (using distance)
input.similarity<-as.data.frame(1-input.ibs$dist)
input.similarity<-as.data.frame(1/(1+input.ibs$dist))

input.similarity<-input.similarity[sample.order,sample.order]

pheatmap(input.similarity, scale = "none", 
         cluster_rows = F, cluster_cols = F,
         display_numbers = T,
         fontsize = 6,
         cellwidth = 16, cellheight = 16,
         width = 4.5, height = 4,
         angle_col = 45,
         filename = paste0("aory",samplegroup,".similarity_matrix.png"))

input.pair<-snpgdsPairScore(input.gds,
                            "AoH","AoL",
                            method = "IBS",
                            type=c("per.pair"),
                            with.id = T)

rm(input.cuttree, input.diss, input.gds, input.ibs, input.similarity, treefile)

##handle vcf file 
samplegroup<-".new"
samplegroup<-".all"
input.vcf<-read.delim(paste0("aory",samplegroup,".vcf"), 
                      header = F, comment.char = "#")

####new
names(input.vcf)<-c("CHROM","POS","ID",
                    "REF","ALT","QUAL",
                    "FILTER","INFO","FORMAT",
                    "AoH","AoL",
                    "RD3230418","RDS170320",
                    "TS3230418","TSS190621",
                    "koji20231118","koji20231209",
                    "koji20230603","koji20230527")

input.vcf<-input.vcf[,c("CHROM","POS","ID",
                        "REF","ALT",
                        "QUAL","FILTER",
                        "INFO","FORMAT",
                        "AoL","TS3230418","TSS190621",
                        "AoH","RD3230418","RDS170320",
                        "koji20230527","koji20230603",
                        "koji20231118","koji20231209")]
####all
names(input.vcf)<-c("CHROM","POS","ID",
                    "REF","ALT","QUAL","FILTER","INFO","FORMAT",
                    "AoH","AoL",
                    "RD3230418","RDS170320",
                    "TS3230418","TSS190621",
                    "koji20231118","koji20231209",
                    "koji2018A","koji2018B",
                    "koji2018C","koji2018D",
                    "koji20230603","koji20230527")

input.vcf<-input.vcf[,c("CHROM","POS","ID",
                      "REF","ALT",
                      "QUAL","FILTER",
                      "INFO","FORMAT",
                      "AoL","TS3230418","TSS190621",
                      "AoH","RD3230418","RDS170320",
                      "koji2018A","koji2018B",
                      "koji2018C","koji2018D",
                      "koji20230527","koji20230603",
                      "koji20231118","koji20231209")]

var.freq<-input.vcf[,c("CHROM","POS","REF","ALT")]
var.summary<-as.data.frame(matrix(NA,
                                  nrow = ncol(input.vcf)-9, 
                                  ncol = 4))
names(var.summary)<-c("sample","SNP","INDEL","unmut")

var.dp<-as.data.frame(matrix(NA,
                             nrow = 0, 
                             ncol = 4))
names(var.dp)<-c("CHROM","POS","REF","ALT")

for (i in 10:ncol(input.vcf)) {
  var.sub<-input.vcf[,c(1,2,4,5,i)]
  var.sub<-separate(var.sub, colnames(input.vcf)[i], 
                    c("GT","AD","DP","GQ","PL"),
                    sep = ":", remove = TRUE, convert = TRUE)
  var.sub<-separate(var.sub, "AD", 
                    c("Dref","Dalt1","Dalt2","Dalt3","Dalt4"),
                    sep = ",", remove = FALSE, convert = TRUE)
  var.sub$Dref[is.na(var.sub$Dref)==T]=20
  var.sub$Dalt1[is.na(var.sub$Dalt1)==T]=0
  var.sub$Dalt2[is.na(var.sub$Dalt2)==T]=0
  var.sub$Dalt3[is.na(var.sub$Dalt3)==T]=0
  var.sub$Dalt4[is.na(var.sub$Dalt4)==T]=0
  var.sub.dp<-var.sub[is.na(var.sub$Dref)==T | (var.sub$Dref+var.sub$Dalt1+var.sub$Dalt2+var.sub$Dalt3+var.sub$Dalt4)>10,]
  var.dp<-rbind(var.dp,var.sub.dp[,c("CHROM","POS","REF","ALT")])
  var.summary$sample[i-9]<-colnames(input.vcf)[i]
  var.summary$SNP[i-9]<-nrow(var.sub.dp[nchar(var.sub.dp$REF)==1 & nchar(var.sub.dp$ALT)==1 & is.na(var.sub.dp$AD)==F,])
  var.summary$INDEL[i-9]<-nrow(var.sub.dp[(nchar(var.sub.dp$REF)!=1 | nchar(var.sub.dp$ALT)!=1) & is.na(var.sub.dp$AD)==F,])
  var.summary$unmut[i-9]<-nrow(var.sub.dp[is.na(var.sub.dp$AD)==T,])
  
  var.sub$freq<-var.sub$Dref/(var.sub$Dref+var.sub$Dalt1+var.sub$Dalt2+var.sub$Dalt3+var.sub$Dalt4)
  var.sub<-var.sub[,c("CHROM","POS","REF","ALT","freq")]
  names(var.sub)[names(var.sub)=="freq"]<-colnames(input.vcf)[i]
  var.freq<-merge(var.freq,var.sub, 
                  by=c("CHROM","POS","REF","ALT"),all=T)
}
rm(var.sub,var.sub.dp,i)
var.freq$CHROM<-as.numeric(gsub("Contig","",var.freq$CHROM))
var.freq<-var.freq[order(var.freq$CHROM,var.freq$POS),]
var.freq$CHROM<-paste0("Contig",var.freq$CHROM)
write.table(var.freq, paste0("aory",samplegroup,".reffq.txt"),
            sep = "\t",
            quote = F, row.names = F)

write.table(var.summary, paste0("aory",samplegroup,".summary.txt"),
            sep = "\t",
            quote = F, row.names = F)

var.dp<-unique(var.dp)
var.dp<-var.dp[order(var.dp$CHROM,var.dp$POS),]
write.table(var.dp, paste0("aory",samplegroup,".vardp10.txt"),
            sep = "\t",
            quote = F, row.names = F)

##pearson correlation
var.freq.matrix<-var.freq
var.freq.matrix$CHROM<-as.numeric(gsub("Contig","",var.freq.matrix$CHROM))
row.names(var.freq.matrix)<-paste0(var.freq.matrix$CHROM,":",var.freq.matrix$POS,":",
                                   var.freq.matrix$REF,"|",var.freq.matrix$ALT)
row.anno<-data.frame(Contig = factor(var.freq.matrix[,c("CHROM")]))
row.names(row.anno)<-row.names(var.freq.matrix)
var.freq.matrix<-var.freq.matrix[,5:ncol(var.freq.matrix)]

d<-cor(var.freq.matrix, method = "pearson", use = "na.or.complete")
d<-as.data.frame(d)
write.table(d, paste0("aory",samplegroup,".pearson.txt"),
            sep = "\t",
            quote = F, row.names = T)

pheatmap(var.freq.matrix, scale = "none", 
         cluster_rows = F, cluster_cols = F,
         display_numbers = F,
         cellwidth = 15, cellheight = 0.14,
         show_rownames = F, show_colnames = T,
         annotation_row = row.anno,
         angle_col = 45,
         width = ncol(var.freq.matrix)*0.3+1, height = nrow(var.freq.matrix)*0.002+1,
         filename = paste0("aory",samplegroup,".varheatmapchr.png"))
pheatmap(var.freq.matrix, scale = "none", 
         cluster_rows = T, cluster_cols = F,
         display_numbers = F,
         cellwidth = 15, cellheight = 0.1,
         show_rownames = F, show_colnames = T,
         angle_col = 45,
         width = ncol(var.freq.matrix)*0.3+1, height = nrow(var.freq.matrix)*0.002+1,
         filename = paste0("aory",samplegroup,".varheatmap.png"))

##filter markers
var.select<-var.freq[var.freq$AoH > 0.95 & var.freq$RD3230418 > 0.95 & var.freq$RDS170320 > 0.95 &
                       var.freq$AoL < 0.05 & var.freq$TS3230418 < 0.05 & var.freq$TSS190621 < 0.05 & 
                       var.freq$koji20230527 > 0.05 & var.freq$koji20230603 > 0.05 &
                       var.freq$koji20231118 > 0.05 & var.freq$koji20231209 > 0.05,]
var.select<-merge(var.select, var.dp, 
                  by = c("CHROM","POS","REF","ALT"), all = F)
write.table(var.select, paste0("aory",samplegroup,".varselect.txt"),
            sep = "\t",
            quote = F, row.names = F)

var.select.vcf<-merge(input.vcf, var.select[,c("CHROM","POS","REF","ALT")],
                      by = c("CHROM","POS","REF","ALT"), all = F)
write.table(var.select.vcf, paste0("aory",samplegroup,".varselectvcf.txt"),
            sep = "\t",
            quote = F, row.names = F)

var.select$CHROM<-as.numeric(gsub("Contig","",var.select$CHROM))
var.select<-var.select[order(var.select$CHROM,var.select$POS),]

var.select.matrix<-var.select
row.names(var.select.matrix)<-paste0(var.select.matrix$CHROM,":",var.select.matrix$POS,":",
                                     var.select.matrix$REF,"|",var.select.matrix$ALT)
row.anno<-data.frame(Contig = factor(var.select.matrix[,c("CHROM")]))
row.names(row.anno)<-row.names(var.select.matrix)

var.select.matrix<-var.select.matrix[,5:ncol(var.select.matrix)]
pheatmap(var.select.matrix, scale = "none", 
         cluster_rows = F, cluster_cols = F,
         display_numbers = F,
         annotation_row = row.anno,
         cellwidth = 50, cellheight = 6,
         fontsize_row = 6, fontsize_col = 12,
         show_rownames = T, show_colnames = T,
         angle_col = 45,
         width = ncol(var.select.matrix)+1, height = nrow(var.select.matrix)*0.1+1,
         filename = paste0("aory",samplegroup,".varselect.heatmap.png"))

##95 percent confidence interval
ci.summary<-as.data.frame(matrix(NA, ncol = 4, nrow = ncol(var.select.matrix)))
names(ci.summary)<-c("Sample","Mean","CI.lower","CI.upper")
ci.summary$Sample<-names(var.select.matrix)
for (i in 1:ncol(var.select.matrix)) {
  t.test.sub<-t.test(var.select.matrix[i])
  ci.summary$Mean[i]<-round(as.numeric(t.test.sub$estimate),3)
  ci.summary$CI.lower[i]<-round(as.numeric(t.test.sub$conf.int)[1],3)
  ci.summary$CI.upper[i]<-round(as.numeric(t.test.sub$conf.int)[2],3)
}

