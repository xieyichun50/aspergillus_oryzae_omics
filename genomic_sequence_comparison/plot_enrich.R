library(ggplot2)
library(dplyr)
library(tidyr)

addline_format <- function(x,...){ 
  gsub('_','\n',x) 
} 

kog2name<-read.delim("/home/yichun/3enrichment/kog2name.txt", header = TRUE)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

ko2name<-read.delim("/home/yichun/3enrichment/ko2name.txt")
kegg2name<-read.delim("/home/yichun/3enrichment/kegg2name.txt")

kegg2ont<- read.delim("/home/yichun/3enrichment/kegglevel.AC.txt", 
                      sep = "\t", colClasses = "character")
names(kegg2ont)[2]="ID"

go2name<-read.delim("/home/yichun/3enrichment/go2name.txt", 
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

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
plotin<-rbind(aoh,aol)
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
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,0.05))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  facet_grid(ONTOLOGY~Cluster, scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a
ggsave(paste0(ref,"KOG.tiff"), width = 12, height = 6, units = "in", dpi = 300)
ggsave(paste0(ref,"KOG.png"), width = 12, height = 6, units = "in", dpi = 300)

##KEGG
anno<-"KEGG"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
plotin<-rbind(aoh,aol)
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

plotdata$Groups<-gsub(" VS ", "\nVS\n", plotdata$Groups)
plotdata1<-subset(plotdata, 
                  #Count > 3 & ratio2>0 & p.adjust< 0.2 & 
                    ONTOLOGY != "Human\nDiseases" & ONTOLOGY != "Organismal\nSystems" & 
                    Description != "MAPK signaling pathway - fly" &
                    Description != "Apoptosis - fly" &
                    Description != "Autophagy - animal" &
                    Description != "Autophagy - other" &
                    Description != "Drug metabolism - other enzymes" &
                    Description != "Mitophagy - animal")
a<-plotdata1 %>%
  #mutate(Groups = factor(Groups, levels = group.order.sub)) %>% 
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12),
        strip.text.x = element_text(size = 11, angle = 0),
        strip.text.y = element_text(size = 9, angle = 0))+
  facet_grid(ONTOLOGY~Cluster, scales = "free", space = "free")

a
ggsave(paste0(ref,"KEGG.tiff"), width = 12, height = 16, units = "in", dpi = 300)
ggsave(paste0(ref,"KEGG.png"), width = 12, height = 16, units = "in", dpi = 300)

##KO
anno="KO"
aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
plotin<-rbind(aoh,aol)
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

plotdata1<-subset(plotdata, ratio1>=0.5)

plotdata1$treat<-paste0(plotdata1$Cluster,plotdata1$Groups)
a<-plotdata1 %>%
  group_by(treat) %>% top_n(-10, p.adjust) %>%
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  facet_grid(~Cluster, scales = "free", space = "free")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12),
        strip.text.x = element_text(size = 11, angle = 0),
        strip.text.y = element_text(size = 9, angle = 0))

a
ggsave(paste0(ref,"ko.tiff"), width = 16, height = 30, units = "in", dpi = 300)
ggsave(paste0(ref,"ko.png"), width = 16, height = 30, units = "in", dpi = 300)

##GO
anno<-"GO"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
plotin<-rbind(aoh,aol)
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"goClass"

plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

GO.L4<-read.delim("/home/yichun/3enrichment/GOlevel4.txt")
plotdata<-plotdata[-which(plotdata$goClass %in% GO.L4$goClass),]

GO.L5<-read.delim("/home/yichun/3enrichment/GOlevel5.txt")
plotdata1<-plotdata[-which(plotdata$goClass %in% GO.L5$goClass),]

GO.L6<-read.delim("/home/yichun/3enrichment/GOlevel6.txt")
plotdata1<-plotdata[which(plotdata$goClass %in% GO.L6$goClass),]

GO.L7<-read.delim("/home/yichun/3enrichment/GOlevel7.txt")
plotdata1<-plotdata[which(plotdata$goClass %in% GO.L7$goClass),]


##all GO
#plotdata1<-subset(plotdata1, ratio1>0.2 & p.adjust < 0.2)
length(unique(plotdata1$goClass))
plotdata1$treat<-paste0(plotdata1$Cluster,plotdata1$Groups)
a<-plotdata1 %>%
  group_by(treat) %>% top_n(-10, p.adjust) %>%
  #mutate(Groups = factor(Groups, levels = group.order.sub)) %>% 
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  facet_grid(ONTOLOGY~Cluster, scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a

ggsave(paste0(ref,"GO.tiff"), width = 12, height = 12, units = "in", dpi = 300)
ggsave(paste0(ref,"GO.png"), width = 12, height = 12, units = "in", dpi = 300)

#top 5 function
plotdata.top5<-plotdata %>% group_by(Cluster) %>% top_n(-10, p.adjust)
#plotdata.top5<-plotdata1 %>% group_by(Cluster) %>% arrange(-ratio2/p.adjust) %>% slice_head(n=5)
write.table(plotdata.top5, 
            paste0(ref,"GO.top5bypadj.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)

##CAZyme
anno<-"CAZyme"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
plotin<-rbind(aoh,aol)
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"

plotdata<-plotinsep

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

list(plotdata$Group)
a<-plotdata[plotdata$kogClass != "3.1.1.72",] %>%
  ggplot(aes(x = Groups, y = kogClass, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  facet_grid(~Cluster, scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a
ggsave(paste0(ref,"CAZyme.tiff"), width = 6, height = 6, units = "in", dpi = 300)
ggsave(paste0(ref,"CAZyme.png"), width = 6, height = 6, units = "in", dpi = 300)

##transporter
anno<-"transporter"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
plotin<-rbind(aoh,aol)
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"

plotdata<-plotinsep

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

list(plotdata$Group)
plotdata<-separate(plotdata, Description, c("Description","ref"), sep = " OS=")
plotdata<-separate(plotdata, Description, c("Description","ref"), sep = " - ")
plotdata<-separate(plotdata, Description, c("Description","ref"), sep = " \\[")
a<-plotdata[plotdata$kogClass != "3.1.1.72",] %>%
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  facet_grid(~Cluster, scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a
ggsave(paste0(ref,"transporter.tiff"), width = 14, height = 6, units = "in", dpi = 300)
ggsave(paste0(ref,"transporter.png"), width = 14, height = 6, units = "in", dpi = 300)


#########AoH REFERENCE
ref<-"Aory.var."

##KOG
anno<-"KOG"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aoh$Groups<-gsub("Unique","AoH.Unique",aoh$Groups)
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
aol<-aol[aol$Groups=="Unique",]
aol$Groups<-gsub("Unique","AoL.Unique",aol$Groups)
plotin<-rbind(aoh,aol)
plotin<-plotin[plotin$Groups %in% c("SNP","INDEL"),]
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"

plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")
plotdata$ONTOLOGY<-gsub(" ","\n",plotdata$ONTOLOGY)

list(plotdata$Group)
a<-plotdata[plotdata$ONTOLOGY != "Poor",] %>%
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,0.05))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12),
        strip.text.y = element_text(size = 9, angle = 0))+
  facet_grid(ONTOLOGY~., scales = "free", space = "free")
a

ggsave(paste0(ref,"KOG.png"), width = 8, height = 4, units = "in", dpi = 300)

##KEGG
anno<-"KEGG"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aoh$Groups<-gsub("Unique","AoH.Unique",aoh$Groups)
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
aol<-aol[aol$Groups=="Unique",]
aol$Groups<-gsub("Unique","AoL.Unique",aol$Groups)
plotin<-rbind(aoh,aol)
plotin<-plotin[plotin$Groups %in% c("SNP","INDEL"),]
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

plotdata$Groups<-gsub(" VS ", "\nVS\n", plotdata$Groups)
plotdata1<-subset(plotdata, 
                  #Count > 3 & ratio2>0 & p.adjust< 0.2 & 
                  ONTOLOGY != "Human\nDiseases" & ONTOLOGY != "Organismal\nSystems" & 
                    Description != "MAPK signaling pathway - fly" &
                    Description != "Apoptosis - fly" &
                    Description != "Autophagy - animal" &
                    Description != "Autophagy - other" &
                    Description != "Drug metabolism - other enzymes" &
                    Description != "Mitophagy - animal")
a<-plotdata1 %>%
  #mutate(Groups = factor(Groups, levels = group.order.sub)) %>% 
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12),
        strip.text.x = element_text(size = 11, angle = 0),
        strip.text.y = element_text(size = 9, angle = 0))+
  facet_grid(ONTOLOGY~., scales = "free", space = "free")

a
ggsave(paste0(ref,"KEGG.png"), width = 8, height = 16, units = "in", dpi = 300)

##KO
anno="KO"
aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aoh$Groups<-gsub("Unique","AoH.Unique",aoh$Groups)
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
aol<-aol[aol$Groups=="Unique",]
aol$Groups<-gsub("Unique","AoL.Unique",aol$Groups)
plotin<-rbind(aoh,aol)
plotin<-plotin[plotin$Groups %in% c("SNP","INDEL"),]
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

plotdata<-subset(plotinsep, is.na(Description) == FALSE & p.adjust < 1)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

plotdata1<-subset(plotdata, ratio1>=0.5)

plotdata1$treat<-paste0(plotdata1$Cluster,plotdata1$Groups)
a<-plotdata1 %>%
  group_by(treat) %>% top_n(-10, p.adjust) %>%
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12),
        strip.text.x = element_text(size = 11, angle = 0),
        strip.text.y = element_text(size = 9, angle = 0))

a

ggsave(paste0(ref,"ko.png"), width = 10, height = 12, units = "in", dpi = 300)

##GO
anno<-"GO"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aoh$Groups<-gsub("Unique","AoH.Unique",aoh$Groups)
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
aol<-aol[aol$Groups=="Unique",]
aol$Groups<-gsub("Unique","AoL.Unique",aol$Groups)
plotin<-rbind(aoh,aol)
plotin<-plotin[plotin$Groups %in% c("SNP","INDEL"),]
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"goClass"

plotdata<-merge(plotinsep, go2name, by = c("goClass"), all.x = TRUE)

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

GO.L4<-read.delim("/home/yichun/3enrichment/GOlevel4.txt")
plotdata<-plotdata[-which(plotdata$goClass %in% GO.L4$goClass),]

GO.L5<-read.delim("/home/yichun/3enrichment/GOlevel5.txt")
plotdata1<-plotdata[-which(plotdata$goClass %in% GO.L5$goClass),]

GO.L6<-read.delim("/home/yichun/3enrichment/GOlevel6.txt")
plotdata1<-plotdata[which(plotdata$goClass %in% GO.L6$goClass),]

GO.L7<-read.delim("/home/yichun/3enrichment/GOlevel7.txt")
plotdata1<-plotdata[which(plotdata$goClass %in% GO.L7$goClass),]


##all GO
#plotdata1<-subset(plotdata1, ratio1>0.2 & p.adjust < 0.2)
length(unique(plotdata1$goClass))
plotdata1$treat<-paste0(plotdata1$Cluster,plotdata1$Groups)
a<-plotdata1 %>%
  group_by(treat) %>% top_n(-10, p.adjust) %>%
  #mutate(Groups = factor(Groups, levels = group.order.sub)) %>% 
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  facet_grid(ONTOLOGY~., scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a

ggsave(paste0(ref,"GO.png"), width = 8, height = 6, units = "in", dpi = 300)

#top 5 function
plotdata.top5<-plotdata %>% group_by(Cluster) %>% top_n(-10, p.adjust)
#plotdata.top5<-plotdata1 %>% group_by(Cluster) %>% arrange(-ratio2/p.adjust) %>% slice_head(n=5)
write.table(plotdata.top5, 
            paste0(ref,"GO.top5bypadj.txt"),
            sep = '\t', row.names = FALSE,
            quote = FALSE)

##CAZyme
anno<-"CAZyme"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aoh$Groups<-gsub("Unique","AoH.Unique",aoh$Groups)
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
aol<-aol[aol$Groups=="Unique",]
aol$Groups<-gsub("Unique","AoL.Unique",aol$Groups)
plotin<-rbind(aoh,aol)
plotin<-plotin[plotin$Groups %in% c("SNP","INDEL"),]
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"

plotdata<-plotinsep

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

list(plotdata$Group)
a<-plotdata[plotdata$kogClass != "3.1.1.72",] %>%
  ggplot(aes(x = Groups, y = kogClass, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  theme(strip.text = element_text(size = 8))
a

ggsave(paste0(ref,"CAZyme.png"), width = 6, height = 4, units = "in", dpi = 300)

##transporter
anno<-"transporter"

aoh<-read.delim(paste0("AoH.",anno,".table.txt"))
aoh$Cluster<-"AoH"
aoh$Groups<-gsub("Unique","AoH.Unique",aoh$Groups)
aol<-read.delim(paste0("AoL.",anno,".table.txt"))
aol$Cluster<-"AoL"
aol<-aol[aol$Groups=="Unique",]
aol$Groups<-gsub("Unique","AoL.Unique",aol$Groups)
plotin<-rbind(aoh,aol)
plotin<-plotin[plotin$Groups %in% c("SNP","INDEL"),]
plotin<-subset(plotin, is.na(Description)==FALSE)

plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)

colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"

plotdata<-plotinsep

plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
plotdata$Group=paste0(plotdata$Cluster, " (", plotdata$Genedenominator, ")")

list(plotdata$Group)
plotdata<-separate(plotdata, Description, c("Description","ref"), sep = " OS=")
plotdata<-separate(plotdata, Description, c("Description","ref"), sep = " - ")
plotdata<-separate(plotdata, Description, c("Description","ref"), sep = " \\[")
a<-plotdata[plotdata$kogClass != "3.1.1.72",] %>%
  ggplot(aes(x = Groups, y = Description, size = Count, colour = ratio1))+
  labs(title = "", size = "Gene count", colour = "Gene ratio", xlab = "", ylab = "")+
  geom_point(shape = 19)+scale_size_area()+
  scale_color_gradient(high = "#FF0000", low = "#0000FF", limits = c(0,1))+
  guides(size = guide_legend(order = 1))+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        axis.title.x = element_text(colour = "white"), 
        axis.title.y = element_text(colour = "white"), 
        axis.text = element_text(colour = "black", size = 11), 
        axis.text.x = element_text(size = 11), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 11), 
        plot.title = element_text(size = 12))+
  #facet_grid(~Cluster, scales = "free", space = "free")+
  theme(strip.text = element_text(size = 8))
a

ggsave(paste0(ref,"transporter.png"), width = 8, height = 4, units = "in", dpi = 300)
