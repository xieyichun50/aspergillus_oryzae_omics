library(dplyr)
library(tidyr)
library(stringr)

setwd("annotation/")
ref<-"AoH"
ref<-"AoL"

##GO
pathways<-read.delim(paste0(ref,".go.anno.xls"), header = T)
pathways<-pathways[,c(1,2)]
names(pathways)<-c("Genes","GO")
GenesGOpair.1v1<-as.data.frame(matrix(NA, nrow = 0, ncol = 2))
names(GenesGOpair.1v1)[1]="Genes"
names(GenesGOpair.1v1)[2]="GO"

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$GO[1], ' G')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(GO, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  GenesGOpair.1v1<-rbind(GenesGOpair.1v1, pairtable.new)
}

GenesGOpair.1v1<-GenesGOpair.1v1[grep("O:", GenesGOpair.1v1$GO),]
GenesGOpair.1v1$GO<-gsub("O:", "GO:", GenesGOpair.1v1$GO)
GenesGOpair.1v1<-separate(GenesGOpair.1v1, GO, 
                          c("GO", "Description", "ONTOLOGY"),
                          sep = ";")

GenesGOpair.1v1<-subset(GenesGOpair.1v1, 
                        is.na(GenesGOpair.1v1$Genes)==FALSE,
                        select = c("GO", "Genes", "Description", "ONTOLOGY"))
GenesGOpair.1v1<-unique(GenesGOpair.1v1)

write.table(GenesGOpair.1v1, 
            file = paste0(ref,".GO.1v1.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)

##KEGG
pathways<-read.delim(paste0(ref, ".kegg.anno.xls"), header = T)
pathways<-pathways[pathways$Ko_class!="",c("Gene_id","Ko_class")]
names(pathways)<-c("Genes","KEGG")

GenesKEGGpair.1v1<-matrix(NA, nrow = 0, ncol = 2)
GenesKEGGpair.1v1<-as.data.frame(GenesKEGGpair.1v1)
names(GenesKEGGpair.1v1)[1]="Genes"
names(GenesKEGGpair.1v1)[2]="KEGG"

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$KEGG[1], '] ')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(KEGG, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  GenesKEGGpair.1v1<-rbind(GenesKEGGpair.1v1, pairtable.new)
}
GenesKEGGpair.1v1<-separate(GenesKEGGpair.1v1, KEGG,
                            c("ONTOLOGY","B","Description"), 
                            sep = ";")
GenesKEGGpair.1v1<-separate(GenesKEGGpair.1v1, Description,
                            c("Description","KEGG"), 
                            sep = " \\[")
GenesKEGGpair.1v1$KEGG<-gsub("PATH:","",GenesKEGGpair.1v1$KEGG)
GenesKEGGpair.1v1$KEGG<-gsub("\\]","",GenesKEGGpair.1v1$KEGG)

GenesKEGGpair.1v1<-subset(GenesKEGGpair.1v1, 
                          is.na(GenesKEGGpair.1v1$Genes)==FALSE, 
                          select = c("KEGG", "Genes",
                                     "Description", "ONTOLOGY","B"))
GenesKEGGpair.1v1<-unique(GenesKEGGpair.1v1)

write.table(GenesKEGGpair.1v1, 
            file = paste0(ref,".KEGG.1v1.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)

##KO
pathways<-read.delim(paste0(ref, ".kegg.anno.xls"), header = T)
pathways<-pathways[pathways$Ko_defi!="",c("Ko_id","Gene_id","Ko_defi")]
names(pathways)<-c("ko","Genes","Description")

write.table(pathways, 
            file = paste0(ref,".KO.1v1.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pathways)

##KOG
pathways<-read.delim(paste0(ref,".kog.anno.xls"), header = T)
pathways<-pathways[pathways$Functional_class!="",c("Gene_id","Functional_class")]
names(pathways)[1]="Genes"
names(pathways)[2]="KOG"

#Format GenesKOGpair
GenesKOGpair.1v1<-matrix(NA, nrow = 0, ncol = 2)
GenesKOGpair.1v1<-as.data.frame(GenesKOGpair.1v1)
names(GenesKOGpair.1v1)[1]="Genes"
names(GenesKOGpair.1v1)[2]="KOG"

for (i in 1:nrow(pathways)) {
  subtable<-pathways[i,]
  rcnames<-list(c(strsplit(subtable$Genes[1], ',')[[1]]),c(strsplit(subtable$KOG[1], '')[[1]]))
  pairtable<-matrix(data = NA, nrow = length(rcnames[[1]]), ncol = length(rcnames[[2]]), dimnames = rcnames)
  pairtable<-as.data.frame(pairtable)
  pairtable$Genes<-rownames(pairtable)
  rownames(pairtable)<-1:nrow(pairtable)
  pairtable<-as.data.frame(pairtable)
  pairtable.new<-pairtable %>% gather(KOG, pair, c(1:ncol(pairtable)-1))
  pairtable.new<-pairtable.new[,c(1:2)]
  GenesKOGpair.1v1<-rbind(GenesKOGpair.1v1, pairtable.new)
}
GenesKOGpair.1v1<-subset(GenesKOGpair.1v1, 
                         is.na(GenesKOGpair.1v1$Genes)==FALSE, 
                         select = c("KOG", "Genes"))
GenesKOGpair.1v1<-unique(GenesKOGpair.1v1)

write.table(GenesKOGpair.1v1, 
            file = paste0(ref,".KOG.1v1.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)
rm(pairtable, pairtable.new, rcnames, subtable, pathways)

##CAZyme
pathways<-read.delim(paste0(ref,".cazy.anno.xls"), header = F, skip = 1)
pathways<-pathways[,c("V6","V1","V7")]
names(pathways)<-c("CAZyme","Genes", "Description")
pathways<-unique(pathways)

write.table(pathways, 
            file = paste0(ref,".CAZyme.1v1.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)

##Transporter
pathways<-read.delim(paste0(ref,".tcdb.anno.xls"), header = T)
pathways<-pathways[,c("TC_number", "Gene_id", "TCDB_gene_discription")]
names(pathways)<-c("transporter","Genes","Description")
pathways<-unique(pathways)

write.table(pathways, 
            file = paste0(ref,".transporter.1v1.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)

