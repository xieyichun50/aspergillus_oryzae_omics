library(dplyr)
library(tidyr)
library(stringr)
setwd("/mnt/content_176/yichun/fungi/aspergillus/comparative/variant")

ref<-"AoH-AoL"
treat<-"indel"
##Indel
snp<-read.delim(paste0(ref,".delta.",treat,".anno.vcf"), header = F, skip = 5)
snp<-separate(snp,V8,c("ANN","Effect"), sep = ";")
snp<-separate(snp,Effect,c("Effect"), sep = ",")
snp<-separate(snp,Effect, 
              c("Allele","Annotation", "Annotation_Impact",
                "Gene_Name","Gene_ID","Feature_Type","Feature_ID",
                "Transcript_BioType","Rank",
                "HGVS.c","HGVS.p",
                "cDNA.pos / cDNA.length","CDS.pos / CDS.length",
                "AA.pos / AA.length","Distance"),
              sep = "\\|")
snp.summary<-as.data.frame(xtabs(~Annotation+Annotation_Impact, snp))
write.table(snp, paste0(ref,".",treat,".effect.xls"),
            row.names = F, quote = F, sep = "\t")
length(unique(snp[snp$Annotation_Impact=="HIGH","Gene_ID"]))
unique(snp[snp$Annotation_Impact=="MODIFIER","Annotation"])

snp.filter<-snp[-which(snp$Annotation %in% c("intergenic_region", 
                                             "upstream_gene_variant",
                                             "downstream_gene_variant",
                                             "intron_variant")),]
snp.filter.summary<-as.data.frame(xtabs(~Annotation, snp.filter))

snp.gene<-as.data.frame(unique(snp.filter[,c("Gene_ID","Annotation","Annotation_Impact")]))
names(snp.gene)[1]<-"Genes"
write.table(snp.gene, paste0(ref,".delta.",treat,".genelist.txt"),
            row.names = F, quote = F, sep = "\t")

##SNP
treat<-"snp"
snp<-read.delim(paste0(ref,".delta.",treat,".anno.vcf"), header = F, skip = 5)
snp<-separate(snp,V8,c("Effect"), sep = ",")
snp<-separate(snp,Effect, 
              c("Allele","Annotation", "Annotation_Impact",
                "Gene_Name","Gene_ID","Feature_Type","Feature_ID",
                "Transcript_BioType","Rank",
                "HGVS.c","HGVS.p",
                "cDNA.pos / cDNA.length","CDS.pos / CDS.length",
                "AA.pos / AA.length","Distance"),
              sep = "\\|")
snp.summary<-as.data.frame(xtabs(~Annotation+Annotation_Impact, snp))

write.table(snp, paste0(ref,".",treat,".effect.xls"),
            row.names = F, quote = F, sep = "\t")

snp.filter<-snp[-which(snp$Annotation %in% c("intergenic_region", 
                                             "upstream_gene_variant",
                                             "downstream_gene_variant",
                                             "intron_variant")),]
length(unique(snp.filter$Gene_ID))
snp.filter.summary<-as.data.frame(xtabs(~Annotation, snp.filter))

snp.gene<-as.data.frame(unique(snp.filter[,c("Gene_ID","Annotation","Annotation_Impact")]))
names(snp.gene)[1]<-"Genes"
write.table(snp.gene, paste0(ref,".delta.",treat,".genelist.txt"),
            row.names = F, quote = F, sep = "\t")
