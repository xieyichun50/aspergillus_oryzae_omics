dir_fq=/mnt/content_176/yichun/fungi/aspergillus/202210AoryzaeRNA/N2207409_80-969971807_RNA/cleandata
dir_new_fq=/mnt/content_176/yichun/fungi/aspergillus/transcriptome/fastq
cat samples_n_reads_decribed.txt | while read spc id1 id2 SRR;
do
echo "fastp -i ${dir_fq}/${SRR}_clean_R1.fq.gz -I ${dir_fq}/${SRR}_clean_R2.fq.gz -o ${dir_new_fq}/$id2.r1.fq -O ${dir_new_fq}/$id2.r2.fq --cut_right -W 4 -M 20 -l 15 -j ${dir_new_fq}/$id2.fastp.json -h ${dir_new_fq}/$id2.fastp.html -R "$id2.fastp_report" -w 38"
fastp -i ${dir_fq}/${SRR}_clean_R1.fq.gz -I ${dir_fq}/${SRR}_clean_R2.fq.gz -o ${dir_new_fq}/$id2.r1.fq -O ${dir_new_fq}/$id2.r2.fq --cut_right -W 4 -M 20 -l 15 -j ${dir_new_fq}/$id2.fastp.json -h ${dir_new_fq}/$id2.fastp.html -R "$id2.fastp_report" -w 38
done
