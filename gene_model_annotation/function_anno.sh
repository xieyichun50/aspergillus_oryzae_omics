/usr/bin/Rscript /mnt/content_176/yichun/scripts/comparative/4function_anno2tree/01gene2function.1v1.R -i AoH.eggnog.emapper.annotations
/usr/bin/Rscript /mnt/content_176/yichun/scripts/comparative/4function_anno2tree/01gene2function.1v1.R -i AoL.eggnog.emapper.annotations
/usr/bin/Rscript /mnt/content_176/yichun/scripts/comparative/4function_anno2tree/02orthogroups2function.R -i /mnt/content_176/yichun/fungi/aspergillus/comparative/OrthoFinder/Orthogroups.tsv

##get the following
all.Genesorthogrouppair.1v1.txt  all.Orthogroups.GO.txt  all.Orthogroups.KEGG.txt  all.Orthogroups.KOG.txt
AoH.Genesorthogrouppair.1v1.txt  AoH.GO.1v1.txt    AoH.ko.1v1.txt    AoH.KEGG.1v1.txt  AoH.KOG.1v1.txt
AoL.Genesorthogrouppair.1v1.txt  AoL.GO.1v1.txt    AoL.ko.1v1.txt    AoL.KEGG.1v1.txt  AoL.KOG.1v1.txt
