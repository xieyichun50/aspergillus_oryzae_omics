library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbreak)
library(eoffice)

kog2name<-read.delim("/mnt/content_176/yichun/scripts/comparative/4function_anno2tree/annotation_generalterm/kog2name.txt",
                     header = TRUE)
kog2name<-kog2name[order(kog2name$ONTOLOGY,kog2name$kogClass),]
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage and processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes and signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

kog.aoh<-read.delim("genome/AoH.KOG.1v1.txt", header = T)
kog.aol<-read.delim("genome/AoL.KOG.1v1.txt", header = T)

kog.aoh.freq<-as.data.frame(xtabs(~KOG, kog.aoh))
kog.aoh.freq$spc<-"AoH"
kog.aol.freq<-as.data.frame(xtabs(~KOG, kog.aol))
kog.aol.freq$spc<-"AoL"

kog.freq<-rbind(kog.aoh.freq, kog.aol.freq)
names(kog.freq)[names(kog.freq)=="KOG"]<-"kogClass"
rm(kog.aoh,kog.aoh.freq,kog.aol,kog.aol.freq)

kog.freq<-merge(kog.freq,kog2name, by = "kogClass", all.x = T)
kog.freq<-kog.freq[kog.freq$kogClass != "S",]
p<- kog.freq %>%
  mutate(spc = factor(spc, levels = c("AoL","AoH"))) %>% 
  ggplot(aes(x = Freq, y = kogName, fill = spc))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), size=3.5, hjust = 0)+
  scale_fill_manual(values = c("grey","black"))+
  scale_x_continuous(limits = c(0,1000))+
  labs(title = "", 
       x = "Gene Count", 
       y = "",
       fill = "")+
  geom_vline(xintercept = 0, color = "black", linetype = "solid", size = 0.5)+
  facet_grid(ONTOLOGY~., scales = "free", space = "free")+
  theme(axis.line = element_line(linetype = "solid"),
        axis.ticks.y = element_line(colour = "black", size = 0.5),
        axis.ticks.x = element_line(colour = "black", size = 0.5),
        axis.text.x = element_text(size = 9, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 11, colour = "black"),
        axis.title.x = element_text(size = 9, colour = "black"),
        panel.background = element_rect(colour = "black", fill = NA),
        panel.grid=element_blank(),
        legend.background = element_rect(fill = NA),
        legend.text = element_text(size = 11, colour = "black"),
        legend.position = c(0.9,0.95))
p
ggsave("figure/KOG.AoH-AoL.png", width = 12, height = 12, units = "in", dpi = 300)
f = "figure/KOG.AoH-AoL.pptx"
topptx(p,f, width = 12, height = 12, units = "in")
