library(dplyr)
library(tidyr)
library(ggbreak)
library(ggplot2)
library(clusterProfiler)

setwd("/mnt/content_176/yichun/fungi/aspergillus/")

##dn/ds statistics
{
genes.NS<-read.delim("genome/AoH_VS_AoL.NSfinal.Comeron.e10.txt", header = T)
names(genes.NS)[names(genes.NS)=="dNdS"]<-"NS"
genes.NS<-separate(genes.NS, "Genes", c("Genes"), sep = "-")
genes.hour<-genes.NS[is.na(genes.NS$NS)==F,c("Genes","NS")]

p<-ggplot(data = genes.hour, aes(x = log2(NS+1), y = stat(density)/sum(stat(density))))+
  geom_histogram(color = "black", fill = "black", binwidth = 0.01)+
  scale_x_continuous(limits = c(0,4), breaks = seq(0,4,1))+
  scale_y_continuous(limits = c(0,0.03),
                     breaks = seq(0,0.03, 0.005))+
  #scale_x_break(c(7, 10))+
  geom_vline(xintercept = 1, color = "red", linetype = "dashed")+
  annotate(geom = "text", 
           label="Negative (Purifying)\nSelection",
           x=0.4, y=0.03,
           colour="black", size=3, fontface = "italic")+
  annotate(geom = "text", 
           label="Positive (Darwinian)\nSelection",
           x=1.6, y=0.03,
           colour="black", size=3, fontface = "italic")+
  labs(x = "log2(dN/dS+1)", y = "Ratio")+
  theme(plot.subtitle = element_text(vjust = 1),
        plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 0, hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.line = element_line(size = 0.5),
        panel.background = element_rect(fill = NA),
        legend.position = "none")

p

ggsave("genome/AoH_VS_AoL.NSfinal.Comeron.e10.png", 
       width = 6, height = 4.5, units = "in", dpi = 300)
}

##dN/dS function annotation --AoH
genes.hour$Groups[genes.hour$NS < 1]<-"Negative selection"
genes.hour$Groups[genes.hour$NS > 1]<-"Positive selection"

kog2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kog2name.txt",
                     header = T)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

GenesKOGpair.1v1<-read.delim("genome/AoH.KOG.1v1.txt", header = T)
genes.hour$Groups[genes.hour$NS < 0.2]<-"Negative selection 5"
genes.hour$Groups[genes.hour$NS >= 0.2 & genes.hour$NS < 0.4]<-"Negative selection 4"
genes.hour$Groups[genes.hour$NS >= 0.4 & genes.hour$NS < 0.6]<-"Negative selection 3"
genes.hour$Groups[genes.hour$NS >= 0.6 & genes.hour$NS < 0.8]<-"Negative selection 2"
genes.hour$Groups[genes.hour$NS >= 0.8 & genes.hour$NS < 1]<-"Negative selection 1"
genes.hour$Groups[genes.hour$NS > 1 & genes.hour$NS < 1.15]<-"Positive selection 1"
genes.hour$Groups[genes.hour$NS >= 1.15 & genes.hour$NS < 1.3]<-"Positive selection 2"
genes.hour$Groups[genes.hour$NS >= 1.3 & genes.hour$NS < 1.5]<-"Positive selection 3"
genes.hour$Groups[genes.hour$NS >= 1.5 & genes.hour$NS < 2]<-"Positive selection 4"
genes.hour$Groups[genes.hour$NS >= 2]<-"Positive selection 5"

KOG.all.1<-compareCluster(Genes ~ Groups, 
                          data = genes.hour, 
                          fun = 'enricher',
                          TERM2GENE = GenesKOGpair.1v1,
                          TERM2NAME = kog2name,
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 1,
                          minGSSize = 1,
                          maxGSSize = 200000)
KOG.all.1<-as.data.frame(KOG.all.1)

write.table(KOG.all.1, "genome/AoH_VS_AoL.NSfinal.Comeron.e10.KOG.table.txt",
            quote = F, row.names = F, sep = "\t")

plotin<-subset(KOG.all.1, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"

plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")
plotdata$Groups<-gsub(" ","\n",plotdata$Groups)
list(plotdata$Group)
a<-plotdata[plotdata$ONTOLOGY != "Poor" & plotdata$ONTOLOGY != "Function unknown",] %>%
  mutate(Groups = factor(Groups, 
                         levels = c("Negative\nselection\n5",
                                    "Negative\nselection\n4",
                                    "Negative\nselection\n3",
                                    "Negative\nselection\n2",
                                    "Negative\nselection\n1",
                                    "Positive\nselection\n1",
                                    "Positive\nselection\n2",
                                    "Positive\nselection\n3",
                                    "Positive\nselection\n4",
                                    "Positive\nselection\n5"))) %>%
  ggplot(aes(x = Groups, y = Description, size = ratio2, colour = ratio1))+
  labs(title = "", size = "Gene Ratio\n \nOver stratum", colour = "Over genome", 
       x = "Purifying <<                                                                                                                                        >> Diverging", y = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,0.2))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.text = element_text(colour = "black", size = 11), 
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  facet_grid(ONTOLOGY~., scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a
ggsave("genome/AoH_VS_AoL.NSfinal.Comeron.e10.KOG.png", width = 16, height = 8, units = "in", dpi = 300)
