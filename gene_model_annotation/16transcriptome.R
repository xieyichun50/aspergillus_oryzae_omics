library(dplyr)
library(tidyr)
library(Hmisc)
library(stringr)
library(ggplot2)
library(eoffice)
library(pheatmap)
library(pathview)
library(clusterProfiler)
setwd("transcriptome/")

species="AoH"
species="AoL"
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
species="single"
filelist<-c("edgeR_gene.min_reps2.min_cpm1_AoH/gene_count_matrix_AoH.csv.AoHES_vs_AoHME.edgeR.DE_results",
            "edgeR_gene.min_reps2.min_cpm1_AoH/gene_count_matrix_AoH.csv.AoHMS_vs_AoHES.edgeR.DE_results",
            "edgeR_gene.min_reps2.min_cpm1_AoL/gene_count_matrix_AoL.csv.AoLES_vs_AoLME.edgeR.DE_results",
            "edgeR_gene.min_reps2.min_cpm1_AoL/gene_count_matrix_AoL.csv.AoLMS_vs_AoLES.edgeR.DE_results")

DEG<-as.data.frame(matrix(NA, nrow = 0, ncol = 7))
names(DEG)<-c("Genes","sampleA","sampleB","logFC","logCPM","PValue","FDR")
for (i in 1:length(filelist)) {
  DEG.sub<-read.delim(file = filelist[i], header = TRUE)
  DEG.sub$Genes<-row.names(DEG.sub)
  row.names(DEG.sub)<-1:nrow(DEG.sub)
  DEG<-rbind(DEG, DEG.sub)
}
rm(DEG.sub,i)

DEG$species<-"AoH"
DEG$species[DEG$sampleA %in% c("AoLME","AoLES","AoLMS")]<-"AoL"
DEG$sampleA<-gsub("AoH","",DEG$sampleA)
DEG$sampleA<-gsub("AoL","",DEG$sampleA)
DEG$sampleB<-gsub("AoH","",DEG$sampleB)
DEG$sampleB<-gsub("AoL","",DEG$sampleB)
DEG$Groups<-paste0(DEG$sampleB,"->",DEG$sampleA)
DEG$change<-"ns"
DEG$change[DEG$FDR < 0.05 & DEG$logFC > 1]<-"Up"
DEG$change[DEG$FDR < 0.05 & DEG$logFC < -1]<-"Down"

DEG.summary<-as.data.frame(xtabs(~species+Groups+change,DEG))
write.table(DEG.summary, 
            paste0("analysis/DEG.summary.",species,".txt"), 
                                row.names = F, quote = F, sep = "\t")
piep<-DEG.summary %>%
  mutate(Groups = factor(Groups, levels = c("ME->ES", "ES->MS")),
         change = factor(change, levels = c("Up","ns","Down"))) %>%
  ggplot(aes(x="", y=Freq, fill=change))+
  geom_bar(stat = "identity", position = position_fill())+
  scale_fill_manual(values=c("coral", "gray", "lightskyblue"))+
  labs(fill = "Change")+
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), 
            color = "black", size = 2)+
  coord_polar(theta = "y", start=0)+
  facet_wrap(species~Groups)+
  theme_void()+
  theme(strip.text = element_text(size = 11))

piep

ggsave(paste0("analysis/DEG.pie.",species,".png"),
       width = 5, height = 5, units = "in", dpi = 300)

####volcano plot
a<-DEG %>%
  mutate(Groups = factor(Groups, levels = c("ME->ES", "ES->MS")),
         change = factor(change, levels = c("Up","ns","Down"))) %>%
  ggplot(aes(x = -logFC, y = -log10(FDR), 
                          colour = change))+
  geom_point(shape = 20, size = 1)+
  scale_colour_manual(values = c("lightskyblue", "gray", "coral"))+
  scale_x_continuous(limits = c(-10,10), breaks = seq(-10,10,2))+
  scale_y_continuous(limits = c(0,20), breaks = seq(0,20,5))+
  labs(x = "log2(Fold change)", y = "-log10(FDR)", legend = "")+
  geom_hline(yintercept = -log10(0.05), color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = -1, color = "grey30", linetype = "dashed")+
  facet_wrap(species~Groups)+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 8))
a
ggsave(paste0("analysis/DEG.Volcano.",species,".png"), width = 5, height = 5, units = "in", dpi = 300)

##heatmap
species="AoH"
species="AoL"
exp.table<-read.delim(file = paste0("gene_count_matrix.TMM_log2_",species,".xls"), 
                      header = TRUE)
names(exp.table)<-gsub("AoH","",names(exp.table))
names(exp.table)<-gsub("AoL","",names(exp.table))
exp.table<-exp.table[,c("ME.rep1", "ME.rep2", "ME.rep3",
                        "ES.rep1", "ES.rep2", "ES.rep3",
                        "MS.rep1", "MS.rep2", "MS.rep3")]
exp.table<-exp.table[row.names(exp.table) %in% unique(DEG$Genes[DEG$species == species & 
                                                                  DEG$change != "ns"]),]

bk.limit<-1
bk <- c(seq(-bk.limit,-0.1,by=0.1),seq(0,bk.limit,by=0.1))
heatp<-pheatmap(exp.table,
                color = c(colorRampPalette(colors = c("lightskyblue","yellow"))(length(bk)/2),
                          colorRampPalette(colors = c("yellow","coral"))(length(bk)/2)),
                legend_breaks=seq(-bk.limit,bk.limit,5), breaks=bk,
                legend = TRUE,
                scale = "row",
                cluster_rows = TRUE, cluster_cols = FALSE,
                cutree_rows = 6, treeheight_row = 30,
                border_color = NA,
                cellwidth = 20, cellheight = 0.075,
                angle_col = c("45"),
                width = 10, height = 10,
                show_colnames = TRUE,
                show_rownames = FALSE,
                filename = paste0("analysis/DEG.heatmap.",species,".png"))
heatp.clust<-cbind(exp.table, 
                   order = as.list(heatp$tree_row["order"]),
                   cluster = cutree(heatp$tree_row, k=6))

##Functional enrich
ref = "AoH"
ref = "AoL"
DEG.sub<-DEG[DEG$species == ref & 
                        DEG$change != "ns",]
##KOG
kog2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kog2name.txt",
                     header = T)
GenesKOGpair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".KOG.1v1.txt"), header = T)
KOG.all.1<-compareCluster(Genes ~ Groups+change, 
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
write.table(KOG.all.1, paste0("analysis/KOG.table.",ref,".txt"),
            quote = F, row.names = F, sep = "\t")

##KEGG
kegg2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kegg2name.txt",
                      header = T)
GenesKEGGpair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".KEGG.1v1.txt"), header = T)
KEGG.all.1<-compareCluster(Genes ~ Groups+change, 
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
write.table(KEGG.all.1, paste0("analysis/KEGG.table.",ref,".txt"),
            quote = F, row.names = F, sep = "\t")

##KO
ko2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/ko2name.txt",
                    header = T)
Geneskopair.1v1<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/",ref,".ko.1v1.txt"), header = T)
ko.all.1<-compareCluster(Genes ~ Groups+change, 
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
write.table(ko.all.1, paste0("analysis/KO.table.",ref,".txt"),
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
GO.all.1<-compareCluster(Genes ~ Groups+change, 
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
write.table(GO.all.1, paste0("analysis/GO.table.",ref,".txt"),
            quote = F, row.names = F, sep = "\t")
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
plotin.sub1<-read.delim(paste0("analysis/KOG.table.","AoH",".txt"), header = T)
plotin.sub1$species="AoH"
plotin.sub2<-read.delim(paste0("analysis/KOG.table.","AoL",".txt"), header = T)
plotin.sub2$species="AoL"
plotin<-rbind(plotin.sub1,plotin.sub2)

plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"
plotinsep$Group=paste0(plotinsep$Groups, ".", plotinsep$species)
plotblank<-as.data.frame(expand_grid(kogClass=unique(plotinsep$kogClass),
                                     Group=unique(plotinsep$Group),
                                     change=unique(plotinsep$change)))
plotdata<-merge(plotblank,plotinsep,by = c("kogClass","Group","change"), all.x = T)
plotdata<-merge(plotdata, kog2name, by = c("kogClass"), all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator

plotdata$ratio1[plotdata$change=="Down"]<-(0-plotdata$ratio1[plotdata$change=="Down"])

plotdata$ratio1[is.na(plotdata$ratio1)==T]<-0

a<-plotdata[plotdata$ONTOLOGY != "Poor",] %>%
  mutate(Group = factor(Group, levels = c("ES->MS.AoL",
                                          "ES->MS.AoH",
                                          "ME->ES.AoL", 
                                          "ME->ES.AoH"))) %>%
  ggplot(aes(x = ratio1, y = kogName, fill = Group))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = c("lightskyblue","dodgerblue4","coral", "coral4"))+
  #scale_x_continuous(limits = c(0,1000))+
  labs(title = "", 
       x = "Down-regulated                Gene ratio                Up-regulated", 
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

ggsave("analysis/KOG.table.png", width = 12, height = 12, units = "in", dpi = 300)

##KEGG
plotin.sub1<-read.delim(paste0("analysis/KEGG.table.","AoH",".txt"), header = T)
plotin.sub1$species="AoH"
plotin.sub2<-read.delim(paste0("analysis/KEGG.table.","AoL",".txt"), header = T)
plotin.sub2$species="AoL"
plotin<-rbind(plotin.sub1,plotin.sub2)

plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

plotinsep$Group=paste0(plotinsep$Groups, ".", plotinsep$species)
plotinsep<-plotinsep[plotinsep$ID %in% plotinsep$ID[plotinsep$p.adjust < 0.2],]
plotblank<-as.data.frame(expand_grid(ID=unique(plotinsep$ID),
                                     Group=unique(plotinsep$Group),
                                     change=unique(plotinsep$change)))
plotblank<-merge(plotblank,unique(plotinsep[,c("ID","Description")]), by = "ID", all.x = T)
plotdata<-merge(plotblank,plotinsep,by = c("ID","Description","Group","change"), all.x = T)
plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator

plotdata$ratio1[plotdata$change=="Down"]<-(0-plotdata$ratio1[plotdata$change=="Down"])

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
  mutate(Group = factor(Group, levels = c("ES->MS.AoL",
                                          "ES->MS.AoH",
                                          "ME->ES.AoL", 
                                          "ME->ES.AoH"))) %>%
  ggplot(aes(x = ratio1, y = Description, fill = Group))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = c("lightskyblue","dodgerblue4","coral", "coral4"))+
  scale_x_continuous(limits = c(-1,1))+
  labs(title = "", 
       x = "Down-regulated                Gene ratio                Up-regulated", 
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
        legend.position = c(0.8,0.8))
a

ggsave("analysis/KEGG.table.png", width = 12, height = 14, units = "in", dpi = 300)

##pathview
Geneskopair.1v1.h<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/AoH.ko.1v1.txt"), header = T)
Geneskopair.1v1.l<-read.delim(paste0("/mnt/content_176/yichun/fungi/aspergillus/genome/AoL.ko.1v1.txt"), header = T)
Geneskopair.1v1<-rbind(Geneskopair.1v1.h,Geneskopair.1v1.l)
rm(Geneskopair.1v1.h,Geneskopair.1v1.l)
DEG.ko.matrix<-as.data.frame(unique(Geneskopair.1v1$ko[Geneskopair.1v1$Genes %in% DEG$Genes]))
names(DEG.ko.matrix)[1]="ko"
DEG$Group<-paste0(DEG$Groups, ".", DEG$species)
Group<-c("ME->ES.AoH","ME->ES.AoL", "ES->MS.AoH","ES->MS.AoL")
DEG.matrix<-merge(DEG,Geneskopair.1v1, by="Genes", all.x = T)

for (i in 1:length(Group)) {
  DEG.matrix.sub<-DEG.matrix[DEG.matrix$Group==Group[i],c("logFC","ko")]
  DEG.ko.matrix.sub<-DEG.matrix.sub %>% group_by(ko) %>% summarise(sum = sum(logFC))
  names(DEG.ko.matrix.sub)[2]=Group[i]
  DEG.ko.matrix<-merge(DEG.ko.matrix, DEG.ko.matrix.sub, all.x = TRUE, by = "ko") 
}
 
rm(DEG.matrix.sub, DEG.ko.matrix.sub)

DEG.ko.matrix<-DEG.ko.matrix[is.na(DEG.ko.matrix$ko)==F,]
row.names(DEG.ko.matrix)<-DEG.ko.matrix$ko
DEG.ko.matrix<-DEG.ko.matrix[,names(DEG.ko.matrix) %in% Group]

write.table(DEG.ko.matrix, "analysis/KO.matrix.txt",
            sep = "\t", quote = F, row.names = T)
setwd("analysis/pv/")
for (i in 93:nrow(keggref)) {
  pathwayid<-keggref$V1[i]
  pathview(gene.data = DEG.ko.matrix,
           pathway.id = pathwayid, 
           species = "ko", 
           gene.idtype = "KEGG", 
           limit = list(gene = 5), 
           bins = list(gene=10), 
           multi.state = TRUE, 
           na.col="transparent", 
           out.suffix = ref)
}

