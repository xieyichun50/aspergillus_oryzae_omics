wd=/mnt/content_176/yichun/fungi/aspergillus/202309aoryzae10
dir_fq=${wd}/fastq
dir_bam=${wd}/alignment_bwa
dir_var=${wd}/variant
refgenome=${dir_bam}/AoH_genomic.fa

samtools faidx ${refgenome}
picard CreateSequenceDictionary --REFERENCE ${refgenome} --OUTPUT ${dir_bam}/AoH_genomic.dict
conda activate gatk

echo "GenomeAnalysisTK -T UnifiedGenotyper -glm BOTH -nt 36 -rf BadCigar --sample_ploidy 1 -R ${refgenome} -o ${dir_var}/aory.merge.vcf -I $(ls ${dir_bam}/*.sorted.bam | sed ':a;N;s/\n/ -I /g;ba;')"
