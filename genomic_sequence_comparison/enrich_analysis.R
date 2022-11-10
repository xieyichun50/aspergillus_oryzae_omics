library(dplyr)
library(clusterProfiler)

ref<-"AoH"
alt<-"AoL"

ref<-"AoL"
alt<-"AoH"
snp<-read.delim(paste0(ref,"-",alt,".delta.snp.genelist.txt"), header = T)
snp$Groups<-"SNP"
indel<-read.delim(paste0(ref,"-",alt,".delta.indel.genelist.txt"), header = T)
indel$Groups<-"INDEL"
unq<-read.delim(paste0(ref,".unique.txt"), header = F)
names(unq)<-"Genes"
unq$Groups<-"Unique"
unq$Genes<-gsub(paste0(ref,"."),"",unq$Genes)

mut<-rbind(snp,indel,unq)

##KOG
kog2name<-read.delim("/Users/XIEYICHUN/Documents/huilab/scripts/eggnog-functional-enrichment-summary-main/general/kog2name.txt",
                     header = T)
GenesKOGpair.1v1<-read.delim(paste0("annotation/",ref,".KOG.1v1.txt"), header = T)
KOG.all.1<-compareCluster(Genes ~ Groups, 
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
kegg2name<-read.delim("/Users/XIEYICHUN/Documents/huilab/scripts/eggnog-functional-enrichment-summary-main/general/kegg2name.txt",
                     header = T)
GenesKEGGpair.1v1<-read.delim(paste0("annotation/",ref,".KEGG.1v1.txt"), header = T)
KEGG.all.1<-compareCluster(Genes ~ Groups, 
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
ko2name<-read.delim("/Users/XIEYICHUN/Documents/huilab/scripts/eggnog-functional-enrichment-summary-main/general/ko2name.txt",
                    header = T)
Geneskopair.1v1<-read.delim(paste0("annotation/",ref,".KO.1v1.txt"), header = T)
ko.all.1<-compareCluster(Genes ~ Groups, 
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
go2name<-read.delim("/Users/XIEYICHUN/Documents/huilab/scripts/eggnog-functional-enrichment-summary-main/general/go2name.txt",
                    header = F)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

GenesGOpair.1v1<-read.delim(paste0("annotation/",ref,".GO.1v1.txt"), header = T)
GO.all.1<-compareCluster(Genes ~ Groups, 
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

##Cazyme
GenesCAZymepair.1v1<-read.delim(paste0("annotation/",ref,".CAZyme.1v1.txt"), header = T)
CAZyme2name<-unique(GenesCAZymepair.1v1[,c("CAZyme","Description")])
CAZyme.all.1<-compareCluster(Genes ~ Groups, 
                         data = mut, 
                         fun = 'enricher',
                         TERM2GENE = GenesCAZymepair.1v1,
                         TERM2NAME = CAZyme2name,
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1,
                         minGSSize = 1,
                         maxGSSize = 2000000)
CAZyme.all.1<-as.data.frame(CAZyme.all.1)
write.table(CAZyme.all.1, paste0(ref,".CAZyme.table.txt"),
            quote = F, row.names = F, sep = "\t")

##Transporter
GenesTranspair.1v1<-read.delim(paste0("annotation/",ref,".transporter.1v1.txt"), header = T)
Transport2name<-unique(GenesTranspair.1v1[,c("transporter","Description")])
Trans.all.1<-compareCluster(Genes ~ Groups, 
                             data = mut, 
                             fun = 'enricher',
                             TERM2GENE = GenesTranspair.1v1,
                             TERM2NAME = Transport2name,
                             pvalueCutoff = 1,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             minGSSize = 1,
                             maxGSSize = 2000000)
Trans.all.1<-as.data.frame(Trans.all.1)
write.table(Trans.all.1, paste0(ref,".transporter.table.txt"),
            quote = F, row.names = F, sep = "\t")
