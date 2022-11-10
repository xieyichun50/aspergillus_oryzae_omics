dir_genome=/mnt/content_176/yichun/fungi/aspergillus/genome
dir_bam=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/bam
dir_gff=/mnt/content_176/yichun/fungi/aspergillus/genome
dir_wd=/mnt/content_176/yichun/fungi/aspergillus/transcriptome

cat samples_n_reads_decribed.txt | cut -f1 | sort | uniq | while read spc;
do
cd ${dir_wd}
cut -f2,3,4 ${dir_wd}/samples_n_reads_decribed_${spc}.txt > ${dir_wd}/samples_n_reads_decribed_${spc}.3.txt
/mnt/content_176/yichun/tools/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix ${dir_wd}/gene_count_matrix_${spc}.csv --method edgeR --samples_file ${dir_wd}/samples_n_reads_decribed_${spc}.3.txt --output ${dir_wd}/edgeR_gene.min_reps2.min_cpm1_${spc} --min_reps_min_cpm 2,1 --contrasts sample_pair_${spc}.txt
cd ${dir_wd}/edgeR_gene.min_reps2.min_cpm1_${spc}
/mnt/content_176/yichun/tools/trinityrnaseq/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ${dir_wd}/gene_count_matrix_${spc}.csv -P 0.05 -C 1 --samples ${dir_wd}/samples_n_reads_decribed_${spc}.txt
done

#${dir_wd}/samples_n_reads_decribed_${spc}.3.txt
AoHES	AoHES.rep1	AoHES1-LEJ5144
AoHES	AoHES.rep2	AoHES2-LEJ5145
AoHES	AoHES.rep3	AoHES4-LEJ5146
AoHME	AoHME.rep1	AoHME3-LEJ5141
AoHME	AoHME.rep2	AoHME4-LEJ5142
AoHME	AoHME.rep3	AoHME5-LEJ5143
AoHMS	AoHMS.rep1	AoHMS1-LEJ5147
AoHMS	AoHMS.rep2	AoHMS2-LEJ5148
AoHMS	AoHMS.rep3	AoHMS4-LEJ5149

#sample_pair_${spc}.txt
AoHES	AoHME
AoHES	AoHMS
AoHME	AoHES
AoHME	AoHMS
AoHMS	AoHES
AoHMS	AoHME
