wd=/mnt/content_176/yichun/fungi/aspergillus/202309aoryzae10
dir_fq=${wd}/fastq
dir_bam=${wd}/alignment
dir_var=${wd}/variant
refgenome=${dir_bam}/AoH_genomic.fa

##call variants
cat ${wd}/samples_n_reads_described.txt | while read id1 id2;
do
bcftools mpileup --threads 36 --max-depth 1000 --min-BQ 20 --min-MQ 20 --excl-flags UNMAP,SECONDARY,QCFAIL,DUP --min-ireads 5 -f ${refgenome} -o ${dir_var}/${id2}.bcf -Ob ${dir_bam}/${id2}.sorted.bam
echo "/mnt/content_176/yichun/fungi/aspergillus/202309aoryzae10/variant/${id2}.bcf" >> ${dir_var}/filelist.txt
done

bcftools merge -l ${dir_var}/filelist.txt --merge both --output-type v --output ${dir_var}/aory.merge.vcf

##create consensus fasta
python3 /mnt/content_176/yichun/tools/vcf2phylip.py --input ${dir_var}/aory.merge.vcf --output-folder ${dir_var}/aory.merge --min-samples-locus 3 --fasta --nexus --write-used-sites
