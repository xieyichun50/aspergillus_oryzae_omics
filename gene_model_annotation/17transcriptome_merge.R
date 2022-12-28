library(dplyr)
library(tidyr)
library(edgeR)
library(Hmisc)
library(stringr)
library(ggplot2)
library(eoffice)
library(pheatmap)
library(pathview)
library(clusterProfiler)
setwd("transcriptome/")

species="AoHL"

####read in expression
exp.table<-read.delim(file = paste0("gene_count_matrix.TMM_log2_",species,".xls"), 
                      header = TRUE)

####read in trait
trait<-read.delim(file = paste0("samples_n_reads_decribed_", species,".txt"), 
                  header = FALSE)

names(trait)=c("species","stage","sample","sampleID")
trait$stage<-gsub("AoH","",trait$stage)
trait$stage<-gsub("AoL","",trait$stage)

##Treatment
exp.mRNA<-exp.table

##PCA
pca.matrix<-prcomp(t(exp.mRNA))
pca.eigenvalues <-pca.matrix$sdev^2
####% of explanation
pca.eigenvalues <- tibble(PC = factor(1:length(pca.eigenvalues)), 
                          variance = pca.eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))

p<-pca.eigenvalues %>% 
  ggplot(aes(x = PC))+
  geom_col(aes(y = pct))+
  geom_line(aes(y = pct_cum, group = 1))+ 
  geom_point(aes(y = pct_cum))+
  scale_y_continuous(breaks = seq(0,100,10))+
  labs(x = "Principal component", y = "Fraction variance explained")
p  
ggsave(paste0("analysis/PCA.variance.sum",species,".png"), width = 8, height = 6, units = "in", dpi = 300)

####distribution
pca.scores <- pca.matrix$x
pca.scores <- pca.scores %>% 
  as_tibble(rownames = "sample")

pca.scores<-merge(pca.scores, trait, by = "sample", all.x = TRUE)
p<-pca.scores %>% 
  mutate(stage = factor(stage, levels = c("ME","ES","MS"))) %>%
  ggplot(aes(x = PC1, y = PC2, shape = factor(stage), fill = factor(species))) +
  geom_point(#alpha = 0.6, 
    size = 3)+
  scale_shape_manual(values = c(21,24,23))+
  scale_fill_manual(values = c("black",
                               "grey"))+
  labs(fill = "isolate", shape = "Stage")+
  guides(fill = guide_legend(order = 1))+
  theme(legend.position = "right")
p
ggsave(paste0("analysis/PCA.PC1-PC2.",species,".png"),
       width = 4, height = 3, units = "in", dpi = 300)


##sample correlation
d<-cor(exp.mRNA, method = "pearson", use = "na.or.complete")
h<-hclust(as.dist(1-d), method = "complete")

png(filename = paste0("analysis/sample_pearson_tree.",species,".png"), 
    res = 300, units = "in", width = 4, height = 6)
plot(varclus(d, similarity="pearson", type = "similarity.matrix"), hang = -1)
dev.off()
rm(pca.eigenvalues,pca.matrix,pca.scores,h,p,d)

## DEG
filelist<-c("edgeR_gene.min_reps2.min_cpm1_AoHL/gene_count_matrix_AoHL.csv.AoLME_vs_AoHME.edgeR.DE_results",
            "edgeR_gene.min_reps2.min_cpm1_AoHL/gene_count_matrix_AoHL.csv.AoLES_vs_AoHES.edgeR.DE_results",
            "edgeR_gene.min_reps2.min_cpm1_AoHL/gene_count_matrix_AoHL.csv.AoLMS_vs_AoHMS.edgeR.DE_results")

DEG<-as.data.frame(matrix(NA, nrow = 0, ncol = 7))
names(DEG)<-c("Genes","sampleA","sampleB","logFC","logCPM","PValue","FDR")
for (i in 1:length(filelist)) {
  DEG.sub<-read.delim(file = filelist[i], header = TRUE)
  DEG.sub$Genes<-row.names(DEG.sub)
  row.names(DEG.sub)<-1:nrow(DEG.sub)
  DEG<-rbind(DEG, DEG.sub)
}
rm(DEG.sub,i)

DEG$stage<-NA
DEG$stage[DEG$sampleA == "AoLME"]<-"ME"
DEG$stage[DEG$sampleA == "AoLES"]<-"ES"
DEG$stage[DEG$sampleA == "AoLMS"]<-"MS"

DEG$sampleA<-gsub("ME","",DEG$sampleA)
DEG$sampleA<-gsub("ES","",DEG$sampleA)
DEG$sampleA<-gsub("MS","",DEG$sampleA)
DEG$sampleB<-gsub("ME","",DEG$sampleB)
DEG$sampleB<-gsub("ES","",DEG$sampleB)
DEG$sampleB<-gsub("MS","",DEG$sampleB)

DEG$Groups<-paste0(DEG$sampleB,"vs",DEG$sampleA)
DEG$change<-"ns"
DEG$change[DEG$FDR < 0.05 & DEG$logFC > 1]<-"AoL"
DEG$change[DEG$FDR < 0.05 & DEG$logFC < -1]<-"AoH"

DEG.summary<-as.data.frame(xtabs(~stage+change,DEG))
write.table(DEG.summary, 
            paste0("analysis/DEG.summary.",species,".txt"), 
            row.names = F, quote = F, sep = "\t")
piep<-DEG.summary %>%
  mutate(Groups = factor(stage, levels = c("ME", "ES", "MS")),
         change = factor(change, levels = c("AoH","ns","AoL"))) %>%
  ggplot(aes(x="", y=Freq, fill=change))+
  geom_bar(stat = "identity", position = position_fill())+
  scale_fill_manual(values=c("coral", "gray", "lightskyblue"))+
  labs(fill = "Change")+
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), 
            color = "black", size = 2)+
  coord_polar(theta = "y", start=0)+
  facet_wrap(stage~.)+
  theme_void()+
  theme(strip.text = element_text(size = 11))

piep

ggsave(paste0("analysis/DEG.pie.",species,".png"),
       width = 8, height = 4, units = "in", dpi = 300)

####volcano plot
a<-DEG %>%
  mutate(stage = factor(stage, levels = c("ME", "ES", "MS")),
         change = factor(change, levels = c("AoH","ns","AoL"))) %>%
  ggplot(aes(x = logFC, y = -log10(FDR), 
             colour = change))+
  geom_point(shape = 20, size = 1)+
  scale_colour_manual(values = c("coral","gray","lightskyblue"))+
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,2))+
  scale_y_continuous(limits = c(0,20), breaks = seq(0,20,5))+
  labs(x = "log2(Fold change)", y = "-log10(FDR)", colour = "High\nexpression")+
  geom_hline(yintercept = -log10(0.05), color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = -1, color = "grey30", linetype = "dashed")+
  facet_wrap(stage~.)+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "right",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 8))
a
ggsave(paste0("analysis/DEG.Volcano.",species,".png"), width = 10, height = 3, units = "in", dpi = 300)

##heatmap
exp.table<-exp.table[,c("AoHME.rep1", "AoHME.rep2", "AoHME.rep3",
                        "AoLME.rep1", "AoLME.rep2", "AoLME.rep3",
                        "AoHES.rep1", "AoHES.rep2", "AoHES.rep3", 
                        "AoLES.rep1", "AoLES.rep2", "AoLES.rep3",
                        "AoHMS.rep1", "AoHMS.rep2", "AoHMS.rep3",
                        "AoLMS.rep1", "AoLMS.rep2", "AoLMS.rep3")]
exp.table<-exp.table[row.names(exp.table) %in% unique(DEG$Genes[DEG$change != "ns"]),]

bk.limit<-1
bk <- c(seq(-bk.limit,-0.1,by=0.1),seq(0,bk.limit,by=0.1))
heatp<-pheatmap(exp.table,
                color = c(colorRampPalette(colors = c("lightskyblue","yellow"))(length(bk)/2),
                          colorRampPalette(colors = c("yellow","coral"))(length(bk)/2)),
                legend_breaks=seq(-bk.limit,bk.limit,5), breaks=bk,
                legend = TRUE,
                scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = 10, treeheight_row = 30,
                border_color = NA,
                cellwidth = 20, cellheight = 0.2,
                angle_col = c("45"),
                width = 10, height = 10,
                show_colnames = TRUE,
                show_rownames = FALSE,
                filename = paste0("analysis/DEG.heatmap.",species,".png"))
heatp.clust<-cbind(exp.table, 
                   order = as.list(heatp$tree_row["order"]),
                   cluster = cutree(heatp$tree_row, k=6))

##Functional enrich
ref="AoH"
DEG.sub<-DEG[DEG$change != "ns",]
##KOG
kog2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kog2name.txt",
                     header = T)
GenesKOGpair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".KOG.1v1.txt"), header = T)
KOG.all.1<-compareCluster(Genes ~ stage+change, 
                          data = DEG.sub, 
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = 200000)
KOG.all.1<-as.data.frame(KOG.all.1)
write.table(KOG.all.1, paste0("analysis/KOG.table.",species,".txt"),
            quote = F, row.names = F, sep = "\t")

##KEGG
kegg2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kegg2name.txt",
                      header = T)
GenesKEGGpair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".KEGG.1v1.txt"), header = T)
KEGG.all.1<-compareCluster(Genes ~ stage+change, 
                           data = DEG.sub, 
                           fun = 'enricher',
                           TERM2GENE = GenesKEGGpair.1v1,
                           TERM2NAME = kegg2name,
                           pvalueCutoff = 1,
                           pAdjustMethod = "BH",
                           qvalueCutoff = 1,
                           minGSSize = 1,
                           maxGSSize = 200000)
KEGG.all.1<-as.data.frame(KEGG.all.1)
write.table(KEGG.all.1, paste0("analysis/KEGG.table.",species,".txt"),
            quote = F, row.names = F, sep = "\t")

##KO
ko2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/ko2name.txt",
                    header = T)
Geneskopair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".ko.1v1.txt"), header = T)
ko.all.1<-compareCluster(Genes ~ stage+change, 
                         data = DEG.sub, 
                         fun = 'enricher',
                         TERM2GENE = Geneskopair.1v1,
                         TERM2NAME = ko2name,
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1,
                         minGSSize = 1,
                         maxGSSize = 200000)
ko.all.1<-as.data.frame(ko.all.1)
write.table(ko.all.1, paste0("analysis/KO.table.",species,".txt"),
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
GO.all.1<-compareCluster(Genes ~ stage+change, 
                         data = DEG.sub, 
                         fun = 'enricher',
                         TERM2GENE = GenesGOpair.1v1,
                         TERM2NAME = go2name,
                         pvalueCutoff = 1,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 1,
                         minGSSize = 1,
                         maxGSSize = 2000000)
GO.all.1<-as.data.frame(GO.all.1)
write.table(GO.all.1, paste0("analysis/GO.table.",species,".txt"),
            quote = F, row.names = F, sep = "\t")

plotin<-GO.all.1
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"goClass"

plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

GO.L4<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/GOlevel4.txt")
plotdata1<-plotdata[-which(plotdata$goClass %in% GO.L4$goClass),]

GO.L5<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/GOlevel5.txt")
plotdata1<-plotdata[-which(plotdata$goClass %in% GO.L5$goClass),]

GO.L6<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/GOlevel6.txt")
plotdata1<-plotdata[which(plotdata$goClass %in% GO.L6$goClass),]

GO.L7<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/GOlevel7.txt")
plotdata1<-plotdata[which(plotdata$goClass %in% GO.L7$goClass),]


##all GO
#top 5 function
plotdata.top5<-plotdata1[plotdata1$p.adjust<0.2,] %>% 
  group_by(Cluster,ONTOLOGY) %>% 
  top_n(-5, ratio1)
write.table(plotdata.top5, 
            paste0("analysis/GO.table.",species,".txt"),
            row.names = F, quote = F,sep = "\t")

##plot

addline_format <- function(x,...){ 
  gsub('_','\n',x) 
} 

kog2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kog2name.txt", 
                     header = TRUE)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

ko2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/ko2name.txt")
kegg2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kegg2name.txt")

kegg2ont<- read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kegglevel.AC.txt", 
                      sep = "\t", colClasses = "character")
names(kegg2ont)[2]="ID"

go2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/go2name.txt", 
                    sep = "\t", colClasses = "character",
                    header = FALSE)
names(go2name)[1]="goClass"
names(go2name)[2]="goName"
names(go2name)[3]="ONTOLOGY"
go2name$ONTOLOGY<-gsub("biological_process", "Biological Process", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("molecular_function", "Molecular Function", go2name$ONTOLOGY)
go2name$ONTOLOGY<-gsub("cellular_component", "Cellular Component", go2name$ONTOLOGY)

##KOG
plotin<-read.delim(paste0("analysis/KOG.table.","AoHL",".txt"), header = T)

plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"
plotinsep$Group=paste0(plotinsep$Groups, ".", plotinsep$species)
plotblank<-as.data.frame(expand_grid(kogClass=unique(plotinsep$kogClass),
                                     stage=unique(plotinsep$stage),
                                     change=unique(plotinsep$change)))
plotdata<-merge(plotblank,plotinsep,by = c("kogClass","stage","change"), all.x = T)
plotdata<-merge(plotdata, kog2name, by = c("kogClass"), all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator

plotdata$ratio1[plotdata$change=="AoH"]<-(0-plotdata$ratio1[plotdata$change=="AoH"])

plotdata$ratio1[is.na(plotdata$ratio1)==T]<-0

a<-plotdata[plotdata$ONTOLOGY != "Poor",] %>%
  mutate(stage = factor(stage, levels = c("ME",
                                          "ES",
                                          "MS"))) %>%
  ggplot(aes(x = ratio1, y = kogName, fill = stage))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = c("lightskyblue","gray","coral"))+
  scale_x_continuous(limits = c(-0.2,0.2))+
  labs(title = "", 
       x = "AoH                Gene ratio                AoL", 
       y = "",
       fill = "")+
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.5)+
  facet_grid(ONTOLOGY~., scales = "free", space = "free")+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11),
        panel.background = element_rect(colour = "black", fill = NA),
        panel.grid=element_blank(),
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 11, colour = "black"),
        legend.position = c(0.88,0.93))
a

ggsave("analysis/KOG.table.AoHL.png", width = 10, height = 8, units = "in", dpi = 300)

##KEGG
plotin<-read.delim(paste0("analysis/KEGG.table.","AoHL",".txt"), header = T)

plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

plotinsep$Group=paste0(plotinsep$Groups, ".", plotinsep$species)
plotinsep<-plotinsep[plotinsep$ID %in% plotinsep$ID[plotinsep$p.adjust < 0.2],]
plotblank<-as.data.frame(expand_grid(ID=unique(plotinsep$ID),
                                     stage=unique(plotinsep$stage),
                                     change=unique(plotinsep$change)))
plotblank<-merge(plotblank,unique(plotinsep[,c("ID","Description")]), by = "ID", all.x = T)
plotdata<-merge(plotblank,plotinsep,by = c("ID","Description","stage","change"), all.x = T)
plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator

plotdata$ratio1[plotdata$change=="AoH"]<-(0-plotdata$ratio1[plotdata$change=="AoH"])

plotdata$ratio1[is.na(plotdata$ratio1)==T]<-0
keggref<-read.delim("/mnt/content_176/yichun/fungi/aspergillus/genome/kegg_aor.txt", header = F)
plotdata1<-plotdata[plotdata$ID %in% keggref$V1 & is.na(plotdata$ONTOLOGY) ==F,]
plotdata1$ONTOLOGY<-gsub("Environmental Information Processing",
                         "Environmental Information\nProcessing",
                         plotdata1$ONTOLOGY)
plotdata1$ONTOLOGY<-gsub("Genetic Information Processing",
                         "Genetic Information\nProcessing",
                         plotdata1$ONTOLOGY)
a<-plotdata1 %>%
  mutate(stage = factor(stage, levels = c("ME",
                                          "ES",
                                          "MS"))) %>%
  ggplot(aes(x = ratio1, y = Description, fill = stage))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = c("lightskyblue","gray","coral"))+
  scale_x_continuous(limits = c(-0.6,0.6))+
  labs(title = "", 
       x = "AoH                Gene ratio                AoL", 
       y = "",
       fill = "")+
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.5)+
  facet_grid(ONTOLOGY~., scales = "free", space = "free")+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 11, colour = "black"),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11),
        strip.text.y = element_text(angle = 0),
        panel.background = element_rect(colour = "black", fill = NA),
        panel.grid=element_blank(),
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 11, colour = "black"),
        legend.position = c(0.8,0.88))
a

ggsave("analysis/KEGG.table.AoHL.png", width = 12, height = 14, units = "in", dpi = 300)

##pathview
Geneskopair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/AoH.ko.1v1.txt"), header = T)
DEG.ko.matrix<-as.data.frame(unique(Geneskopair.1v1$ko[Geneskopair.1v1$Genes %in% DEG$Genes]))
names(DEG.ko.matrix)[1]="ko"
Group<-c("ME","ES", "MS")
DEG.matrix<-merge(DEG,Geneskopair.1v1, by="Genes", all.x = T)

for (i in 1:length(Group)) {
  DEG.matrix.sub<-DEG.matrix[DEG.matrix$stage==Group[i],c("logFC","ko")]
  DEG.ko.matrix.sub<-DEG.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(logFC))
  names(DEG.ko.matrix.sub)[2]=Group[i]
  DEG.ko.matrix<-merge(DEG.ko.matrix, DEG.ko.matrix.sub, all.x = TRUE, by = "ko") 
}

rm(DEG.matrix.sub, DEG.ko.matrix.sub)

DEG.ko.matrix<-DEG.ko.matrix[is.na(DEG.ko.matrix$ko)==F,]
row.names(DEG.ko.matrix)<-DEG.ko.matrix$ko
DEG.ko.matrix<-DEG.ko.matrix[,names(DEG.ko.matrix) %in% Group]

write.table(DEG.ko.matrix, "analysis/KO.matrix.AoHL.txt",
            sep = "\t", quote = F, row.names = T)
setwd("analysis/pv/")
for (i in c(1:90,93:nrow(keggref))) {
  pathwayid<-keggref$V1[i]
  pathview(gene.data = DEG.ko.matrix,
           pathway.id = pathwayid, 
           species = "ko", 
           gene.idtype = "KEGG", 
           limit = list(gene = 5), 
           bins = list(gene=10), 
           multi.state = TRUE, 
           na.col="transparent", 
           out.suffix = species)
}

##match DEG with genetic variants
ref<-"AoH"
alt<-"AoL"

mut<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/comparative/variant/",
                       ref,"-",alt,".delta.mut.genelist.txt"), header = T)
degmut<-merge(DEG.sub, mut, by = "Genes", all = F)
degmut$Annotation.simple<-degmut$Annotation
degmut$Annotation.simple<-gsub("inframe_deletion","synonymous_variant",degmut$Annotation.simple)
degmut$Annotation.simple<-gsub("frameshift_variant","missense_variant",degmut$Annotation.simple)

mut.ko.matrix<-as.data.frame(unique(Geneskopair.1v1$ko[Geneskopair.1v1$Genes %in% degmut$Genes]))
names(mut.ko.matrix)[1]="ko"
mut.ko.matrix.bk<-mut.ko.matrix
mut.matrix<-merge(degmut,Geneskopair.1v1, by = "Genes", all.x = T)
Group<-c("5_prime_UTR_variant",
         "synonymous_variant",
         "missense_variant",
         "gain/loss_start/stop",
         "splicing variant",
         "3_prime_UTR_variant")

for (i in 1:length(Group)) {
  mut.matrix.sub<-mut.matrix[mut.matrix$Annotation==Group[i],
                             c("logFC","ko")]
  mut.ko.matrix.sub<-mut.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(logFC))
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
kegg2ko<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kegg2koID.txt",
                    header = T)
keggref<-merge(mut.ko.matrix,kegg2ko, by = "ko", all.x = T)
keggref<-keggref[,c("KEGG")]
mut.ko.matrix<-mut.ko.matrix[,names(mut.ko.matrix) %in% Group]
for (i in 1:length(keggref)) {
  pathwayid<-keggref[i]
  plotdata<-pathview(gene.data = mut.ko.matrix,
                     pathway.id = pathwayid, 
                     species = "ko", 
                     gene.idtype = "KEGG", 
                     limit = list(gene = 5), 
                     bins = list(gene=10), 
                     multi.state = TRUE, 
                     na.col="transparent", 
                     out.suffix = paste0(ref,".degmut"))
  if (class(plotdata) == "list"){
    plot.data.gene<-as.data.frame(plotdata$plot.data.gene)
    plot.data.gene$kegg.names<-unlist(plot.data.gene$kegg.names)
    plot.data.gene$all.mapped<-unlist(plot.data.gene$all.mapped)
    plot.data.gene$type<-unlist(plot.data.gene$type)
    write.table(plot.data.gene,
                paste0(pathwayid, ".table.degmut.txt"),
                sep = "\t", row.names = F, quote = F)
  } else {}
}

