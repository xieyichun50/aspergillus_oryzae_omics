dir_new_fq=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq
dir_genome=/mnt/content_176/yichun/fungi/aspergillus/genome
dir_bam=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/bam

hisat2-build -p 39 -f ${dir_genome}/AoH_genomic.fa ${dir_genome}/AoH_genomic.fa
hisat2-build -p 39 -f ${dir_genome}/AoL_genomic.fa ${dir_genome}/AoL_genomic.fa
hisat2-build -p 39 -f ${dir_genome}/AoT_genomic.fa ${dir_genome}/AoT_genomic.fa

cat samples_n_reads_decribed.txt | while read spc id1 id2 SRR;
do
echo "hisat2 -x ${dir_genome}/${spc}_genomic.fa -1 ${dir_new_fq}/${id2}.r1.fq -2 ${dir_new_fq}/${id2}.r2.fq -S ${dir_bam}/${id2}.sam --phred33 --dta-cufflinks --novel-splicesite-outfile ${dir_bam}/${id2}.splicesite.txt -p 38 --fr"
hisat2 -x ${dir_genome}/${spc}_genomic.fa -1 ${dir_new_fq}/${id2}.r1.fq -2 ${dir_new_fq}/${id2}.r2.fq -S ${dir_bam}/${id2}.sam --phred33 --dta-cufflinks --novel-splicesite-outfile ${dir_bam}/${id2}.splicesite.txt -p 38 --fr

samtools view --threads 38 -b -S ${dir_bam}/${id2}.sam > ${dir_bam}/${id2}.bam
rm ${dir_bam}/${id2}.sam
samtools sort --threads 38 ${dir_bam}/${id2}.bam -o ${dir_bam}/${id2}.sorted.bam
rm ${dir_bam}/${id2}.bam
samtools index ${dir_bam}/${id2}.sorted.bam
#java -Xmx40g -jar /tools/gatk/picard.jar MarkDuplicates INPUT=$id.sorted.bam OUTPUT=$id.sorted.redu.bam METRICS_FILE=$id.sorted.redu.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
#samtools index $id.sorted.redu.bam
done

