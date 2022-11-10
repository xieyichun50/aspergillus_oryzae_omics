spc=
dir_genome=/mnt/content_176/yichun/fungi/aspergillus/genome
dir_new_fq=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq
dir_bam=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/bam

#funannotate clean -i /mnt/content_176/yichun/fungi/aspergillus/genome/${spc}_genomic.fa -o ${spc}_genomic.clean.fa

echo "funannotate mask --cpu 38 -i /mnt/content_176/yichun/fungi/aspergillus/genome/${spc}_genomic.fa -o ${spc}/${spc}_genomic.masked.fa"

funannotate mask --cpu 38 -i /mnt/content_176/yichun/fungi/aspergillus/genome/${spc}_genomic.fa -o ${spc}/${spc}_genomic.masked.fa

echo "funannotate train -i ${spc}/${spc}_genomic.masked.fa -o ${spc} --left $(find /mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq | xargs ls | grep "/" | grep ${spc} | grep "r1" | tr "\n" " ") --right $(find /mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq | xargs ls | grep "/" | grep ${spc} | grep "r2" | tr "\n" " ") --stranded FR --jaccard_clip --no_trimmomatic --memory 300G --coverage 100 --cpus 76 --species "Aspergillus oryzae" --strain ${spc}"

funannotate train -i ${spc}/${spc}_genomic.masked.fa -o ${spc} --left $(find /mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq | xargs ls | grep "/" | grep ${spc} | grep "r1" | tr "\n" " ") --right $(find /mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq | xargs ls | grep "/" | grep ${spc} | grep "r2" | tr "\n" " ") --stranded FR --jaccard_clip --no_trimmomatic --memory 300G --coverage 100 --cpus 76 --species "Aspergillus oryzae" --strain ${spc}
