##AoH
AoH<-read.delim("/mnt/content_93/home/yichun_hml/aspergillus/rawdata/TGS/03.Genome_Component/AoH/Gene/AoH.gff",
                header = F)
AoH.new<-AoH[-which(AoH$V1 %in%  c("Contig1","Contig2")),]
AoH.adjust<-AoH[AoH$V1 %in%  c("Contig1","Contig2"),]
AoH.adjust$contig<-NA
AoH.adjust$start<-NA
AoH.adjust$end<-NA

AoH.adjust$contig[AoH.adjust$V1=="Contig1" & 
                    AoH.adjust$V4 >0 & 
                    AoH.adjust$V5 <= 3087526]<-"Contig1.1"
AoH.adjust$start[AoH.adjust$V1=="Contig1" & 
                   AoH.adjust$V4 >0 & 
                   AoH.adjust$V5 <= 3087526]<-AoH.adjust$V4[AoH.adjust$V1=="Contig1" & 
                                                              AoH.adjust$V4 >0 & 
                                                              AoH.adjust$V5 <= 3087526]
AoH.adjust$end[AoH.adjust$V1=="Contig1" & 
                 AoH.adjust$V4 >0 & 
                 AoH.adjust$V5 <= 3087526]<-AoH.adjust$V5[AoH.adjust$V1=="Contig1" & 
                                                            AoH.adjust$V4 >0 & 
                                                            AoH.adjust$V5 <= 3087526]

AoH.adjust$contig[AoH.adjust$V1=="Contig1" & 
                    AoH.adjust$V4 >3087526 & 
                    AoH.adjust$V5 <= 5798651]<-"Contig1.2"
AoH.adjust$start[AoH.adjust$V1=="Contig1" & 
                   AoH.adjust$V4 >3087526 & 
                   AoH.adjust$V5 <= 5798651]<-AoH.adjust$V4[AoH.adjust$V1=="Contig1" & 
                                                              AoH.adjust$V4 >3087526 & 
                                                              AoH.adjust$V5 <= 5798651]-3087526
AoH.adjust$end[AoH.adjust$V1=="Contig1" & 
                 AoH.adjust$V4 >3087526 & 
                 AoH.adjust$V5 <= 5798651]<-AoH.adjust$V5[AoH.adjust$V1=="Contig1" & 
                                                            AoH.adjust$V4 >3087526 & 
                                                            AoH.adjust$V5 <= 5798651]-3087526

AoH.adjust$contig[AoH.adjust$V1=="Contig2" & 
                    AoH.adjust$V4 >0 & 
                    AoH.adjust$V5 <= 2406880]<-"Contig2.1"
AoH.adjust$start[AoH.adjust$V1=="Contig2" & 
                   AoH.adjust$V4 >0 & 
                   AoH.adjust$V5 <= 2406880]<-AoH.adjust$V4[AoH.adjust$V1=="Contig2" & 
                                                              AoH.adjust$V4 >0 & 
                                                              AoH.adjust$V5 <= 2406880]
AoH.adjust$end[AoH.adjust$V1=="Contig2" & 
                 AoH.adjust$V4 >0 & 
                 AoH.adjust$V5 <= 2406880]<-AoH.adjust$V5[AoH.adjust$V1=="Contig2" & 
                                                            AoH.adjust$V4 >0 & 
                                                            AoH.adjust$V5 <= 2406880]

AoH.adjust$contig[AoH.adjust$V1=="Contig2" & 
                    AoH.adjust$V4 >2406880 & 
                    AoH.adjust$V5 <= 4430350]<-"Contig2.2"
AoH.adjust$start[AoH.adjust$V1=="Contig2" & 
                   AoH.adjust$V4 >2406880 & 
                   AoH.adjust$V5 <= 4430350]<-AoH.adjust$V4[AoH.adjust$V1=="Contig2" & 
                                                              AoH.adjust$V4 >2406880 & 
                                                              AoH.adjust$V5 <= 4430350]-2406880
AoH.adjust$end[AoH.adjust$V1=="Contig2" & 
                 AoH.adjust$V4 >2406880 & 
                 AoH.adjust$V5 <= 4430350]<-AoH.adjust$V5[AoH.adjust$V1=="Contig2" & 
                                                            AoH.adjust$V4 >2406880 & 
                                                            AoH.adjust$V5 <= 4430350]-2406880

AoH.adjust<-AoH.adjust[,c("contig", "V2","V3",
                          "start", "end",
                          "V6","V7","V8","V9")]
names(AoH.adjust)<-c("V1","V2","V3","V4","V5","V6","V7","V8","V9")

AoH.new<-rbind(AoH.new,AoH.adjust)

write.table(AoH.new, 
            "/mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/update/AoH/AoH.gff",
            quote = F, row.names = F, col.names = F, sep = "\t")

##AoL
AoL<-read.delim("/mnt/content_93/home/yichun_hml/aspergillus/rawdata/TGS/03.Genome_Component/AoL/Gene/AoL.gff",
                header = F)
AoL.new<-AoL[-which(AoL$V1 %in%  c("Contig1")),]
AoL.adjust<-AoL[AoL$V1 %in%  c("Contig1"),]
AoL.adjust$contig<-NA
AoL.adjust$start<-NA
AoL.adjust$end<-NA

AoL.adjust$contig[AoL.adjust$V1=="Contig1" & 
                    AoL.adjust$V4 >0 & 
                    AoL.adjust$V5 <= 868750]<-"Contig1.1"
AoL.adjust$start[AoL.adjust$V1=="Contig1" & 
                   AoL.adjust$V4 >0 & 
                   AoL.adjust$V5 <= 868750]<-AoL.adjust$V4[AoL.adjust$V1=="Contig1" & 
                                                              AoL.adjust$V4 >0 & 
                                                              AoL.adjust$V5 <= 868750]
AoL.adjust$end[AoL.adjust$V1=="Contig1" & 
                 AoL.adjust$V4 >0 & 
                 AoL.adjust$V5 <= 868750]<-AoL.adjust$V5[AoL.adjust$V1=="Contig1" & 
                                                            AoL.adjust$V4 >0 & 
                                                            AoL.adjust$V5 <= 868750]

AoL.adjust$contig[AoL.adjust$V1=="Contig1" & 
                    AoL.adjust$V4 >868750 & 
                    AoL.adjust$V5 <= 3075458]<-"Contig1.2"
AoL.adjust$start[AoL.adjust$V1=="Contig1" & 
                   AoL.adjust$V4 >868750 & 
                   AoL.adjust$V5 <= 3075458]<-AoL.adjust$V4[AoL.adjust$V1=="Contig1" & 
                                                              AoL.adjust$V4 >868750 & 
                                                              AoL.adjust$V5 <= 3075458]-868750
AoL.adjust$end[AoL.adjust$V1=="Contig1" & 
                 AoL.adjust$V4 >868750 & 
                 AoL.adjust$V5 <= 3075458]<-AoL.adjust$V5[AoL.adjust$V1=="Contig1" & 
                                                            AoL.adjust$V4 >868750 & 
                                                            AoL.adjust$V5 <= 3075458]-868750

AoL.adjust$contig[AoL.adjust$V1=="Contig1" & 
                    AoL.adjust$V4 >3075452 & 
                    AoL.adjust$V5 <= 4371685]<-"Contig1.3"
AoL.adjust$start[AoL.adjust$V1=="Contig1" & 
                   AoL.adjust$V4 >3075452 & 
                   AoL.adjust$V5 <= 4371685]<-AoL.adjust$V4[AoL.adjust$V1=="Contig1" & 
                                                              AoL.adjust$V4 >3075452 & 
                                                              AoL.adjust$V5 <= 4371685]-3075452
AoL.adjust$end[AoL.adjust$V1=="Contig1" & 
                 AoL.adjust$V4 >3075452 & 
                 AoL.adjust$V5 <= 4371685]<-AoL.adjust$V5[AoL.adjust$V1=="Contig1" & 
                                                            AoL.adjust$V4 >3075452 & 
                                                            AoL.adjust$V5 <= 4371685]-3075452

AoL.adjust<-AoL.adjust[,c("contig", "V2","V3",
                          "start", "end",
                          "V6","V7","V8","V9")]
names(AoL.adjust)<-c("V1","V2","V3","V4","V5","V6","V7","V8","V9")

AoL.new<-rbind(AoL.new,AoL.adjust)

write.table(AoL.new, 
            "/mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/update/AoL/AoL.gff",
            quote = F, row.names = F, col.names = F, sep = "\t")

