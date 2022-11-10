cat filelist | while read input_fasta input_fastq prefix;
do
minimap2 -a -x map-pb -MD -t 15 ${input_fasta} ${input_fastq} > ${prefix}.sam
echo "Alignment finished!"
samtools view --thread 15 -b -S ${prefix}.sam > ${prefix}.bam
echo "sam to bam conversion finished!"
samtools sort --thread 15 ${prefix}.bam -o ${prefix}.sorted.bam
echo "bam sorted!"
samtools index ${prefix}.sorted.bam
echo "bam indexed!"

samtools flagstat --threads 15 ${prefix}.sorted.bam
echo "Calling SV"
sniffles -m ${prefix}.sorted.bam -v SV.${prefix}.vcf --min_support 10 -t 15

done

echo "Calling SNP"
cat filelist | while read input_fasta input_fastq prefix;
do
echo "bcftools mpileup --thread 15 -Q 7 -f ${input_fasta} -o SNP.${prefix}.vcf -Ov ${prefix}.sorted.bam &"
done 
echo "wait"

##filelist
/mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.fasta	/mnt/content_93/home/yichun_hml/aspergillus/fastq/TGS.AoH.fastq.gz	TGS.AoH.TGS.AoH
/mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.fasta	/mnt/content_93/home/yichun_hml/aspergillus/fastq/TGS.AoL.fastq.gz	TGS.AoH.TGS.AoL
/mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.fasta	/mnt/content_93/home/yichun_hml/aspergillus/fastq/TGS.AoH.fastq.gz	TGS.AoL.TGS.AoH
/mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.fasta	/mnt/content_93/home/yichun_hml/aspergillus/fastq/TGS.AoL.fastq.gz	TGS.AoL.TGS.AoL
