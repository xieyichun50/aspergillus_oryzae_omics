cat sample_info.txt | while read prefix reffa bamfile;
do

#echo "Alignment with pbmm2"
#echo "pbmm2 align --sort -j 40 -J 40 --preset SUBREAD ${reffa} ${bamfile} ${prefix}.aligned.bam"
#pbmm2 align --sort -j 40 -J 40 --preset SUBREAD ${reffa} ${bamfile} ${prefix}.aligned.bam

#pbindex -j 40 ${prefix}.aligned.bam
samtools faidx ${reffa}

echo "/opt/pacbio/smrtlink/smrtcmds/bin/ipdSummary --reference ${reffa} --gff ${prefix}.modified.gff --csv ${prefix}.modified.csv --bigwig ${prefix}.modified.bigwig -j 40 --pvalue 0.01 --identify m6A,m4C --minCoverage 10 --identifyMinCov 15 --methylMinCov 30 -v --methylFraction --outfile ${prefix}.modified ${prefix}.aligned.bam"
/opt/pacbio/smrtlink/smrtcmds/bin/ipdSummary --reference ${reffa} --gff ${prefix}.modified.gff --csv ${prefix}.modified.csv --bigwig ${prefix}.modified.bigwig -j 40 --pvalue 0.01 --identify m6A,m4C --minCoverage 10 --identifyMinCov 15 --methylMinCov 30 -v --methylFraction --outfile ${prefix}.modified ${prefix}.aligned.bam
done
