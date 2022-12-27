library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/mnt/content_176/yichun/fungi/aspergillus/comparative/variant/")
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

ref<-"Aory"

##KOG
anno<-"KOG"

plotin<-read.delim(paste0("AoH.",anno,".table.txt"))
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"

plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

list(plotdata$Group)
a<-plotdata[plotdata$ONTOLOGY != "Poor",] %>%
  mutate(Groups = factor(Groups, levels = c("SNP","INDEL")),
         Annotation = factor(Annotation, levels = c("5_prime_UTR_variant",
                                                    "synonymous_variant",
                                                    "missense_variant",
                                                    "inframe_deletion",
                                                    "frameshift_variant",
                                                    "gain/loss_start/stop",
                                                    "splicing variant",
                                                    "3_prime_UTR_variant"))) %>% 
  ggplot(aes(x = Annotation, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,0.05))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a
ggsave(paste0(ref,".KOG.tiff"), width = 12, height = 8, units = "in", dpi = 300)
ggsave(paste0(ref,".KOG.png"), width = 12, height = 8, units = "in", dpi = 300)

##KEGG
anno<-"KEGG"

plotin<-read.delim(paste0("AoH.",anno,".table.txt"))
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)

plotdata<-merge(plotdata, kegg2ont, by = "ID", all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

plotdata$ONTOLOGY<-gsub(" ", "\n", plotdata$ONTOLOGY)
list(plotdata$Group)

keggref<-read.delim("/mnt/content_176/yichun/fungi/aspergillus/genome/kegg_aor.txt", header = F)
plotdata1<-plotdata[plotdata$ID %in% keggref$V1 & is.na(plotdata$ONTOLOGY) ==F,]
a<-plotdata1[plotdata1$ONTOLOGY != "Human\nDiseases",] %>%
  mutate(Groups = factor(Groups, levels = c("SNP","INDEL")),
         Annotation = factor(Annotation, levels = c("5_prime_UTR_variant",
                                                    "synonymous_variant",
                                                    "missense_variant",
                                                    "inframe_deletion",
                                                    "frameshift_variant",
                                                    "gain/loss_start/stop",
                                                    "splicing variant",
                                                    "3_prime_UTR_variant"))) %>% 
  ggplot(aes(x = Annotation, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12),
        strip.text.x = element_text(size = 11, angle = 0),
        strip.text.y = element_text(size = 9, angle = 0))+
  facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")

a
ggsave(paste0(ref,".KEGG.tiff"), width = 12, height = 16, units = "in", dpi = 300)
ggsave(paste0(ref,".KEGG.png"), width = 12, height = 16, units = "in", dpi = 300)

a<-plotdata1[plotdata1$ONTOLOGY == "Metabolism",] %>%
  mutate(Groups = factor(Groups, levels = c("SNP","INDEL")),
         Annotation = factor(Annotation, levels = c("5_prime_UTR_variant",
                                                    "synonymous_variant",
                                                    "missense_variant",
                                                    "inframe_deletion",
                                                    "frameshift_variant",
                                                    "gain/loss_start/stop",
                                                    "splicing variant",
                                                    "3_prime_UTR_variant"))) %>% 
  ggplot(aes(x = Annotation, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12),
        strip.text.x = element_text(size = 11, angle = 0),
        strip.text.y = element_text(size = 9, angle = 0))+
  facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")

a
ggsave(paste0(ref,".KEGG2.tiff"), width = 12, height = 16, units = "in", dpi = 300)
ggsave(paste0(ref,".KEGG2.png"), width = 12, height = 16, units = "in", dpi = 300)

##KO
anno="KO"

plotin<-read.delim(paste0("AoH.",anno,".table.txt"))
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

plotdata1<-subset(plotdata, ratio1>=0.5)

a<-plotdata1 %>%
  mutate(Groups = factor(Groups, levels = c("SNP","INDEL")),
         Annotation = factor(Annotation, levels = c("5_prime_UTR_variant",
                                                    "synonymous_variant",
                                                    "missense_variant",
                                                    "inframe_deletion",
                                                    "frameshift_variant",
                                                    "gain/loss_start/stop",
                                                    "splicing variant",
                                                    "3_prime_UTR_variant"))) %>% 
  group_by(Cluster) %>% top_n(-5, p.adjust) %>%
  ggplot(aes(x = Annotation, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  facet_grid(~Groups, scales = "free", space = "free")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12),
        strip.text.x = element_text(size = 11, angle = 0),
        strip.text.y = element_text(size = 9, angle = 0))

a
ggsave(paste0(ref,".KO.tiff"), width = 16, height = 30, units = "in", dpi = 300)
ggsave(paste0(ref,".KO.png"), width = 16, height = 30, units = "in", dpi = 300)

##GO
anno<-"GO"

plotin<-read.delim(paste0("AoH.",anno,".table.txt"))
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
plotdata.top5<-plotdata %>% group_by(Cluster) %>% top_n(-10, p.adjust)
#plotdata.top5<-plotdata1 %>% group_by(Cluster) %>% arrange(-ratio2/p.adjust) %>% slice_head(n=5)
write.table(plotdata.top5, 
            paste0(ref,".GO.top5bypadj.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)

plotdata1<-subset(plotdata1, ratio1>0.2 & p.adjust < 0.2)
length(unique(plotdata1$goClass))
a<-plotdata1 %>%
  group_by(Cluster) %>% top_n(-10, p.adjust) %>%
  mutate(Groups = factor(Groups, levels = c("SNP","INDEL")),
         Annotation = factor(Annotation, levels = c("5_prime_UTR_variant",
                                                    "synonymous_variant",
                                                    "missense_variant",
                                                    "inframe_deletion",
                                                    "frameshift_variant",
                                                    "gain/loss_start/stop",
                                                    "splicing variant",
                                                    "3_prime_UTR_variant"))) %>% 
  ggplot(aes(x = Annotation, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  facet_grid(ONTOLOGY~Groups, scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a

ggsave(paste0(ref,".GO.tiff"), width = 18, height = 24, units = "in", dpi = 300)
ggsave(paste0(ref,".GO.png"), width = 18, height = 24, units = "in", dpi = 300)

