/usr/bin/Rscript dnds_calculation.R -q AoH.longest-gene.cds.fa -r AoL.longest-gene.cds.fa -o AoH_VS_AoL -b /home/yichun_hml/miniconda3/bin/blastp -m /home/yichun_hml/miniconda3/bin/mafft | tee AoH_VS_AoL.log

#dnds_calculation.R
##usage
#Rscript dNdS_calculation -q query.cds-transcripts.fa -r reference.cds-transcripts.fa -b /root/miniconda3/bin/blastp -m /root/miniconda3/bin/mafft

#Output the following
##Species.taxidN.txt (unique performed to reduce record numbers)

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(orthologr))

##Create parameters
option_list <- list(
  make_option(c("-q","--query"), type="character", default=NULL,
              help="cds sequence of your targeted species' [default %default]",
              dest="qry"),
  make_option(c("-r","--reference"), type="character", default=NULL,
              help="cds sequence of your reference species' [default %default]",
              dest="ref"),
  make_option(c("-o","--output"), type="character", default=NULL,
              help="output prefix' [default %default]",
              dest="output"),
  make_option(c("-b","--blast"), type="character", default=NULL,
              help="which blastp' [default %default]",
              dest="blast"),
  make_option(c("-m","--mafft"), type="character", default=NULL,
              help="which mafft' [default %default]",
              dest="mafft"),
  make_option(c("-k","--kaks"), type="character", default=NULL,
              help="which KaKs_Calculator' [default %default]",
              dest="kaks"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i orthologs.txt -o out [options]",option_list=option_list)
opt = parse_args(parser)

##manual add files
{
  #opt$qry<-"Coprinopsis_cinerea.longest-gene.cds.fa"
  #opt$ref<-"Fusarium_graminearum.longest-gene.cds.fa"
  #opt$blast<-"/root/miniconda3/bin/blastp"
  #opt$mafft<-"/root/miniconda3/bin/mafft"
  #opt$kaks<-"/root/miniconda3/bin/KaKs_Calculator"
}

#opt$output<-paste0(opt$output,
#                   ".NSfinal.txt")
cat(paste0("File will be written to ",opt$output,"\n"))

# using the `aa_aln_path` or `blast_path` arguments

#rbh.result<-blast_rec(query_file = opt$qry,
#                     subject_file = opt$ref,
#                     seq_type = "cds",
#                     blast_algorithm = "blastp",
#                     eval = "1E-10",
#                     max.target.seqs = 10,
#                     comp_cores = 76,
#                     path = opt$blast,
#                     clean_folders = FALSE,
#                     delete_corrupt_cds = TRUE)
#write.table(rbh.result, paste0(opt$output,".RBH.txt"),
#            row.names = F, sep = "\t", quote = F)

dnds.result<-dNdS(query_file = opt$qry,
                  subject_file = opt$ref,
                  ortho_detection = "RBH",
                  blast_path = opt$blast,
                  eval = "1E-10",
                  aa_aln_type = "multiple",
                  aa_aln_tool = "mafft",
                  aa_aln_path = opt$mafft,
                  codon_aln_tool = "pal2nal",
                  dnds_est.method = "Comeron",
                  kaks_calc_path = opt$kaks,
                  comp_cores = 76,
                  clean_folders = FALSE,
                  delete_corrupt_cds = TRUE)

if (opt$qry=="Coprinopsis_cinerea.longest-gene.cds.fa") {
  IDmatch<-read.delim("Coprinopsis_cinerea.GenematchID", header = F)
  names(IDmatch)<-c("Genes","query_id")
  dnds.result<-merge(dnds.result, IDmatch, by = "query_id", all.x =T)
} else if (opt$qry=="Fusarium_graminearum.longest-gene.cds.fa"){
  IDmatch<-read.delim("Fusarium_graminearum.GenematchID", header = F)
  names(IDmatch)<-c("query_id","Genes")
  dnds.result<-merge(dnds.result, IDmatch, by = "query_id", all.x =T)
} else if (opt$qry=="Neurospora_crassa.longest-gene.cds.fa"){
  IDmatch<-read.delim("Neurospora_crassa.GenematchID", header = F)
  names(IDmatch)<-c("query_id","Genes")
  dnds.result<-merge(dnds.result, IDmatch, by = "query_id", all.x =T)
} else if (opt$qry=="Rhizopus_delemar.longest-gene.cds.fa"){
  IDmatch<-read.delim("Rhizopus_delemar.GenematchID", header = F)
  names(IDmatch)<-c("query_id","Genes")
  dnds.result<-merge(dnds.result, IDmatch, by = "query_id", all.x =T)
} else if (file.exists(paste0(opt$qry,".GenematchID"))){
  IDmatch<-read.delim(paste0(opt$qry,".GenematchID"), header = F)
  names(IDmatch)<-c("query_id","Gene")
  dnds.result<-merge(dnds.result, IDmatch, by = "query_id", all.x =T)
} else {
  dnds.result$Genes<-dnds.result$query_id
}
write.table(dnds.result, paste0(opt$output,".NSfinal.txt"),
            row.names = F, sep = "\t", quote = F)

save.image(paste0(opt$output,".RData"))
