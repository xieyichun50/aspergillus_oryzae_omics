library(dplyr)
library(clusterProfiler)

ref<-"AoH"
alt<-"AoL"

snp<-read.delim(paste0(ref,"-",alt,".delta.snp.genelist.txt"), header = T)
snp$Annotation[grep("gain",snp$Annotation)]<-"gain/loss_start/stop"
snp$Annotation[grep("lost",snp$Annotation)]<-"gain/loss_start/stop"
snp$Annotation[-which(snp$Annotation %in% c("3_prime_UTR_variant",
                                            "5_prime_UTR_variant",
                                            "missense_variant",
                                            "synonymous_variant",
                                            "gain/loss_start/stop"))]<-"splicing variant"
snp$Groups<-"SNP"

indel<-read.delim(paste0(ref,"-",alt,".delta.indel.genelist.txt"), header = T)
indel$Annotation[indel$Annotation %in% c("disruptive_inframe_deletion",
                                         "conservative_inframe_deletion")]<-"inframe_deletion"
indel$Annotation[grep("gain",indel$Annotation)]<-"gain/loss_start/stop"
indel$Annotation[grep("lost",indel$Annotation)]<-"gain/loss_start/stop"
indel$Annotation[-which(indel$Annotation %in% c("3_prime_UTR_variant",
                                                "5_prime_UTR_variant",
                                                "inframe_deletion",
                                                "frameshift_variant",
                                                "gain/loss_start/stop"))]<-"splicing variant"

indel$Groups<-"INDEL"
mut<-rbind(snp,indel)
write.table(mut,
            paste0(ref,"-",alt,".delta.mut.genelist.txt"),
            row.names = F, sep = "\t", quote = F)

##KOG
kog2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kog2name.txt",
                     header = T)
GenesKOGpair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".KOG.1v1.txt"), header = T)
KOG.all.1<-compareCluster(Genes ~ Groups+Annotation, 
                          data = mut, 
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = 200000)
KOG.all.1<-as.data.frame(KOG.all.1)
write.table(KOG.all.1, paste0(ref,".KOG.table.txt"),
            quote = F, row.names = F, sep = "\t")

##KEGG
kegg2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kegg2name.txt",
                      header = T)
GenesKEGGpair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".KEGG.1v1.txt"), header = T)
KEGG.all.1<-compareCluster(Genes ~ Groups+Annotation, 
                           data = mut, 
                           fun = 'enricher',
                           TERM2GENE = GenesKEGGpair.1v1,
                           TERM2NAME = kegg2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 200000)
KEGG.all.1<-as.data.frame(KEGG.all.1)
write.table(KEGG.all.1, paste0(ref,".KEGG.table.txt"),
            quote = F, row.names = F, sep = "\t")

##KO
ko2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/ko2name.txt",
                    header = T)
Geneskopair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".ko.1v1.txt"), header = T)
ko.all.1<-compareCluster(Genes ~ Groups+Annotation, 
                         data = mut, 
                         fun = 'enricher',
                         TERM2GENE = Geneskopair.1v1,
                         TERM2NAME = ko2name,
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1,
                         minGSSize = 1,
                         maxGSSize = 200000)
ko.all.1<-as.data.frame(ko.all.1)
write.table(ko.all.1, paste0(ref,".KO.table.txt"),
            quote = F, row.names = F, sep = "\t")

##GO
go2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/go2name.txt",
                    header = F)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

GenesGOpair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".GO.1v1.txt"), header = T)
GO.all.1<-compareCluster(Genes ~ Groups+Annotation, 
                         data = mut, 
                         fun = 'enricher',
                         TERM2GENE = GenesGOpair.1v1,
                         TERM2NAME = go2name,
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1,
                         minGSSize = 1,
                         maxGSSize = 2000000)
GO.all.1<-as.data.frame(GO.all.1)
write.table(GO.all.1, paste0(ref,".GO.table.txt"),
            quote = F, row.names = F, sep = "\t")

##pathview
keggref<-read.delim("/mnt/content_176/yichun/fungi/aspergillus/genome/kegg_aor.txt", header = F)
Geneskopair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/AoH.ko.1v1.txt"), header = T)
mut.ko.matrix<-as.data.frame(unique(Geneskopair.1v1$ko[Geneskopair.1v1$Genes %in% mut$Genes]))
names(mut.ko.matrix)[1]="ko"
mut.ko.matrix.bk<-mut.ko.matrix
mut$count<-1
mut.matrix<-merge(mut,Geneskopair.1v1, by = "Genes", all.x = T)

Group<-c("SNP","INDEL")

for (i in 1:length(Group)) {
  mut.matrix.sub<-mut.matrix[mut.matrix$Group==Group[i],c("count","ko")]
  mut.ko.matrix.sub<-mut.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(count))
  names(mut.ko.matrix.sub)[2]=Group[i]
  mut.ko.matrix<-merge(mut.ko.matrix, mut.ko.matrix.sub, all.x = TRUE, by = "ko") 
}

rm(mut.matrix.sub, mut.ko.matrix.sub)
mut.ko.matrix<-mut.ko.matrix[is.na(mut.ko.matrix$ko)==F,]
row.names(mut.ko.matrix)<-mut.ko.matrix$ko
mut.ko.matrix<-mut.ko.matrix[,names(mut.ko.matrix) %in% Group]

write.table(mut.ko.matrix, "KO.matrix.txt",
            sep = "\t", quote = F, row.names = T)
mut.ko.matrix.all<-mut.ko.matrix

##SNP
mut.ko.matrix<-mut.ko.matrix.bk
Group<-c("5_prime_UTR_variant",
         "synonymous_variant",
         "missense_variant",
         "gain/loss_start/stop",
         "splicing variant",
         "3_prime_UTR_variant")

for (i in 1:length(Group)) {
  mut.matrix.sub<-mut.matrix[mut.matrix$Annotation==Group[i] & mut.matrix$Groups == "SNP",
                             c("count","ko")]
  mut.ko.matrix.sub<-mut.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(count))
  names(mut.ko.matrix.sub)[2]=Group[i]
  mut.ko.matrix<-merge(mut.ko.matrix, mut.ko.matrix.sub, all.x = TRUE, by = "ko") 
}

rm(mut.matrix.sub, mut.ko.matrix.sub)
mut.ko.matrix<-mut.ko.matrix[is.na(mut.ko.matrix[,2])==F |
                               is.na(mut.ko.matrix[,3])==F |
                               is.na(mut.ko.matrix[,4])==F |
                               is.na(mut.ko.matrix[,5])==F |
                               is.na(mut.ko.matrix[,6])==F |
                               is.na(mut.ko.matrix[,7])==F,]
row.names(mut.ko.matrix)<-mut.ko.matrix$ko
mut.ko.matrix<-mut.ko.matrix[,names(mut.ko.matrix) %in% Group]

write.table(mut.ko.matrix, "KO.matrix.SNP.txt",
            sep = "\t", quote = F, row.names = T)
mut.ko.matrix.SNP<-mut.ko.matrix

##INDEL
mut.ko.matrix<-mut.ko.matrix.bk
Group<-c("5_prime_UTR_variant",
         "inframe_deletion",
         "frameshift_variant",
         "gain/loss_start/stop",
         "splicing variant",
         "3_prime_UTR_variant")
  
for (i in 1:length(Group)) {
  mut.matrix.sub<-mut.matrix[mut.matrix$Annotation==Group[i] & mut.matrix$Groups == "INDEL",
                             c("count","ko")]
  mut.ko.matrix.sub<-mut.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(count))
  names(mut.ko.matrix.sub)[2]=Group[i]
  mut.ko.matrix<-merge(mut.ko.matrix, mut.ko.matrix.sub, all.x = TRUE, by = "ko") 
}

rm(mut.matrix.sub, mut.ko.matrix.sub)

mut.ko.matrix<-mut.ko.matrix[is.na(mut.ko.matrix[,2])==F |
                               is.na(mut.ko.matrix[,3])==F |
                               is.na(mut.ko.matrix[,4])==F |
                               is.na(mut.ko.matrix[,5])==F |
                               is.na(mut.ko.matrix[,6])==F |
                               is.na(mut.ko.matrix[,7])==F,]
row.names(mut.ko.matrix)<-mut.ko.matrix$ko
mut.ko.matrix<-mut.ko.matrix[,names(mut.ko.matrix) %in% Group]

write.table(mut.ko.matrix, "KO.matrix.INDEL.txt",
            sep = "\t", quote = F, row.names = T)
mut.ko.matrix.indel<-mut.ko.matrix

setwd("pv/")
library(pathview)
for (i in c(1:90,93:nrow(keggref))) {
  pathwayid<-keggref$V1[i]
  plotdata<-pathview(gene.data = mut.ko.matrix.all,
                     pathway.id = pathwayid, 
                     species = "ko", 
                     gene.idtype = "KEGG", 
                     limit = list(gene = 5), 
                     bins = list(gene=10), 
                     multi.state = TRUE, 
                     na.col="transparent", 
                     out.suffix = paste0(ref,".all"))
  if (class(plotdata) == "list"){
    plot.data.gene<-as.data.frame(plotdata$plot.data.gene)
    plot.data.gene$kegg.names<-unlist(plot.data.gene$kegg.names)
    plot.data.gene$all.mapped<-unlist(plot.data.gene$all.mapped)
    plot.data.gene$type<-unlist(plot.data.gene$type)
    write.table(plot.data.gene,
                paste0(pathwayid, ".table.txt"),
                sep = "\t", row.names = F, quote = F)
  } else {}
  #SNP
  plotdata<-pathview(gene.data = mut.ko.matrix.SNP,
                     pathway.id = pathwayid, 
                     species = "ko", 
                     gene.idtype = "KEGG", 
                     limit = list(gene = 5), 
                     bins = list(gene=10), 
                     multi.state = TRUE, 
                     na.col="transparent", 
                     out.suffix = paste0(ref,".SNP"))
  if (class(plotdata) == "list"){
    plot.data.gene<-as.data.frame(plotdata$plot.data.gene)
    plot.data.gene$kegg.names<-unlist(plot.data.gene$kegg.names)
    plot.data.gene$all.mapped<-unlist(plot.data.gene$all.mapped)
    plot.data.gene$type<-unlist(plot.data.gene$type)
    write.table(plot.data.gene,
                paste0(pathwayid, ".table.SNP.txt"),
                sep = "\t", row.names = F, quote = F)
  } else {}
  
  plotdata<-pathview(gene.data = mut.ko.matrix.indel,
                     pathway.id = pathwayid, 
                     species = "ko", 
                     gene.idtype = "KEGG", 
                     limit = list(gene = 5), 
                     bins = list(gene=10), 
                     multi.state = TRUE, 
                     na.col="transparent", 
                     out.suffix = paste0(ref,".INDEL"))
  if (class(plotdata) == "list"){
    plot.data.gene<-as.data.frame(plotdata$plot.data.gene)
    plot.data.gene$kegg.names<-unlist(plot.data.gene$kegg.names)
    plot.data.gene$all.mapped<-unlist(plot.data.gene$all.mapped)
    plot.data.gene$type<-unlist(plot.data.gene$type)
    write.table(plot.data.gene,
                paste0(pathwayid, ".table.INDEL.txt"),
                sep = "\t", row.names = F, quote = F)
  } else {}
}
