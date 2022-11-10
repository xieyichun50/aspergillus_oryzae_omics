dir_genome=/mnt/content_176/yichun/fungi/aspergillus/genome
dir_bam=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/bam
dir_gff=/mnt/content_176/yichun/fungi/aspergillus/genome

####Expression level
cat samples_n_reads_decribed.txt | while read spc id1 id2 SRR;
do  
stringtie -e -B -p 38 -G ${dir_gff}/Aspergillus_oryzae_${spc}.gff3 -A ${id2}_gene_count.xls -o ballgown/${spc}/${id2}_ballgown/${id2}.gtf ${dir_bam}/$id2.sorted.bam
done

cat samples_n_reads_decribed.txt | cut -f1 | sort | uniq | while read spc;
do
prepDE.py -i ballgown/${spc} -l 150 -g gene_count_matrix_${spc}.csv -t transcript_count_matrix_${spc}.csv
sed -i 's/,/\t/g;s/_ballgown//g' gene_count_matrix_${spc}.csv transcript_count_matrix_${spc}.csv
Rscript /mnt/content_176/yichun/fungi/hourglass/scripts/edgeR_TMM_norm.R --suffix "_${spc}"
done
