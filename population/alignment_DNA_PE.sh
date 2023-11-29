wd=/mnt/content_176/yichun/fungi/aspergillus/202309aoryzae10
dir_fq=${wd}/fastq
dir_bam=${wd}/alignment
dir_var=${wd}/variant
refgenome=${dir_bam}/AoH_genomic.fa

bowtie2-build -f ${refgenome} ${refgenome}

tail -n2 ${wd}/samples_n_reads_described.txt | while read id1 id2;
do
echo "fastp -i ${dir_fq}/${id2}/*_1.fq* -I ${dir_fq}/${id2}/*_2.fq* -o ${dir_fq}/${id2}.r1.fq -O ${dir_fq}/${id2}.r2.fq --cut_right -W 4 -M 20 -l 15 -j ${dir_fq}/${id2}.fastp.json -h ${dir_fq}/${id2}.fastp.html -R ${id2}.fastp_report -w 16"
fastp -i ${dir_fq}/${id2}/*_1.fq* -I ${dir_fq}/${id2}/*_2.fq* -o ${dir_fq}/${id2}.r1.fq -O ${dir_fq}/${id2}.r2.fq --cut_right -W 4 -M 20 -l 15 -j ${dir_fq}/${id2}.fastp.json -h ${dir_fq}/${id2}.fastp.html -R ${id2}.fastp_report -w 16

echo "bowtie2 -x ${refgenome} -1 ${dir_fq}/${id2}.r1.fq -2 ${dir_fq}/${id2}.r2.fq -S ${dir_bam}/${id2}.sam --phred33 -p 36"
bowtie2 -x ${refgenome} -1 ${dir_fq}/${id2}.r1.fq -2 ${dir_fq}/${id2}.r2.fq -S ${dir_bam}/${id2}.sam --phred33 -p 36

echo "samtools view --threads 36 -b -S ${dir_bam}/${id2}.sam > ${dir_bam}/${id2}.bam"
samtools view --threads 36 -b -S ${dir_bam}/${id2}.sam > ${dir_bam}/${id2}.bam
rm ${dir_bam}/${id2}.sam

echo "samtools sort --threads 36 ${dir_bam}/${id2}.bam -o ${dir_bam}/${id2}.sorted.bam"
samtools sort --threads 36 ${dir_bam}/${id2}.bam -o ${dir_bam}/${id2}.sort.bam
rm ${dir_bam}/${id2}.bam

echo "picard MarkDuplicates -INPUT ${dir_bam}/${id2}.sort.bam -OUTPUT ${dir_bam}/${id2}.sorted.bam -METRICS_FILE ${dir_bam}/${id2}.sort.matrics.txt -REMOVE_DUPLICATES false -ASSUME_SORTED true"
picard MarkDuplicates -INPUT ${dir_bam}/${id2}.sort.bam -OUTPUT ${dir_bam}/${id2}.sorted.bam -METRICS_FILE ${dir_bam}/${id2}.sort.matrics.txt -REMOVE_DUPLICATES false -ASSUME_SORTED true

echo "samtools index -@ 36 ${dir_bam}/${id2}.sorted.bam"
samtools index -@ 36 ${dir_bam}/${id2}.sorted.bam

#bcftools mpileup --threads 36 --max-depth 1000 --min-BQ 20 --min-MQ 20 --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --min-ireads 5 -f ${refgenome} -o ${dir_var}/${id2}.vcf -Ov ${dir_bam}/${id2}.sorted.bam
#bcftools view --types snps  --output-type v --output ${dir_var}/${id2}.snp.vcf ${dir_var}/${id2}.vcf
#bcftools view --types indels --output-type v --output ${dir_var}/${id2}.indel.vcf ${dir_var}/${id2}.vcf

done
