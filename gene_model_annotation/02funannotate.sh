spc=AoH
dir_genome=/mnt/content_176/yichun/fungi/aspergillus/genome
dir_new_fq=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq
dir_bam=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/bam

#echo "funannotate predict -i ${spc}/${spc}_genomic.masked.fa -o ${spc} --species "Aspergillus oryzae" --strain ${spc} --busco_db ascomycota --genemark_mode ET --busco_seed_species "aspergillus_oryzae" --augustus_species "aspergillus_oryzae" --optimize_augustus --cpus 76 --name ${spc}"

#funannotate predict -i ${spc}/${spc}_genomic.masked.fa -o ${spc} --species "Aspergillus oryzae" --strain ${spc} --busco_db ascomycota --genemark_mode ET --busco_seed_species "aspergillus_oryzae" --augustus_species "aspergillus_oryzae" --optimize_augustus --cpus 76 --name ${spc}

#echo "funannotate update -i ${spc} -o ${spc} --jaccard_clip --cpus 76 "

#funannotate update -i ${spc} -o ${spc} --jaccard_clip --cpus 76 

#echo "funannotate iprscan -i ${spc} -m local --cpus 36"
#funannotate iprscan -i ${spc} -m local --cpus 36
#echo "funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m phobius"
#funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m phobius 
echo "funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m antismash"
funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m antismash

#echo "funannotate annotate -i ${spc} --sbt /mnt/content_176/yichun/scripts/genomes/template.sbt --busco_db ascomycota"

#funannotate annotate -i ${spc} --sbt /mnt/content_176/yichun/scripts/genomes/template.sbt --busco_db ascomycota

spc=AoL
dir_genome=/mnt/content_176/yichun/fungi/aspergillus/genome
dir_new_fq=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq
dir_bam=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/bam

#echo "funannotate predict -i ${spc}/${spc}_genomic.masked.fa -o ${spc} --species "Aspergillus oryzae" --strain ${spc} --busco_db ascomycota --genemark_mode ET --busco_seed_species "aspergillus_oryzae" --augustus_species "aspergillus_oryzae" --optimize_augustus --cpus 76 --name ${spc}"

#funannotate predict -i ${spc}/${spc}_genomic.masked.fa -o ${spc} --species "Aspergillus oryzae" --strain ${spc} --busco_db ascomycota --genemark_mode ET --busco_seed_species "aspergillus_oryzae" --augustus_species "aspergillus_oryzae" --optimize_augustus --cpus 76 --name ${spc}

#echo "funannotate update -i ${spc} -o ${spc}  --jaccard_clip --cpus 36"

#funannotate update -i ${spc} -o ${spc} --jaccard_clip --cpus 36

#echo "funannotate iprscan -i ${spc} -m local --cpus 36"
#funannotate iprscan -i ${spc} -m local --cpus 36
#echo "funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m phobius"
#funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m phobius 
echo "funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m antismash"
funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m antismash

#echo "funannotate annotate -i ${spc} --sbt /mnt/content_176/yichun/scripts/genomes/template.sbt --busco_db ascomycota"

#funannotate annotate -i ${spc} --sbt /mnt/content_176/yichun/scripts/genomes/template.sbt --busco_db ascomycota

spc=AoT
dir_genome=/mnt/content_176/yichun/fungi/aspergillus/genome
dir_new_fq=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq
dir_bam=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/bam

#echo "funannotate predict -i ${spc}/${spc}_genomic.masked.fa -o ${spc} --species "Aspergillus oryzae" --strain ${spc} --busco_db ascomycota --genemark_mode ET --busco_seed_species "aspergillus_oryzae" --augustus_species "aspergillus_oryzae" --optimize_augustus --cpus 76 --name ${spc}"

#funannotate predict -i ${spc}/${spc}_genomic.masked.fa -o ${spc} --species "Aspergillus oryzae" --strain ${spc} --busco_db ascomycota --genemark_mode ET --busco_seed_species "aspergillus_oryzae" --augustus_species "aspergillus_oryzae" --optimize_augustus --cpus 76 --name ${spc}

#echo "funannotate update -i ${spc} -o ${spc} --jaccard_clip --cpus 36"

#funannotate update -i ${spc} -o ${spc} --jaccard_clip --cpus 36

#echo "funannotate iprscan -i ${spc} -m local --cpus 36"
#funannotate iprscan -i ${spc} -m local --cpus 36
echo "funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m antismash"
funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m antismash
echo "funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m phobius"
funannotate remote -e xieyichun50@link.cuhk.edu.hk -i ${spc} -m phobius 

#echo "funannotate annotate -i ${spc} --sbt /mnt/content_176/yichun/scripts/genomes/template.sbt --busco_db ascomycota"

#funannotate annotate -i ${spc} --sbt /mnt/content_176/yichun/scripts/genomes/template.sbt --busco_db ascomycota
