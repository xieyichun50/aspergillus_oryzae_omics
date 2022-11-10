library(ape)
library(ggplot2)
library(tidytree)
library(ggtree)
library(flextable)
library(tidyr)
library(dplyr)
library(stringr)
#library(svglite)
library(ggplotify)
library(eoffice)
setwd("phylo/")

##cafe_tree_number
tree<-read.tree("Afla.nwk")
tree.new<-root(tree, "Aspergillus_flavus_NRRL3357.genomic",resolve.root = TRUE)

tree.new<-read.tree("phylo.tree")
p<-ggtree(tree.new, size = 1, layout = "rectangular")+
  geom_tiplab(aes(label=paste("     ",str_replace_all(label,"_"," ")), sep = ""),
              offset = 0.05,
              #fontface = "italic",
              size = 3,
              align=TRUE, linesize=.5)+
  xlim(NA,0.01)

p

ggsave(paste0("phylo-new.png"),
       width = 30, height = 15, units = "in", dpi = 300)

f = paste0("phylo.pptx")
topptx(p,f, width = 15, height = 15, units = "in")
