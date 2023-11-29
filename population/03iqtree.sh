wd=/mnt/content_176/yichun/fungi/aspergillus/202309aoryzae10
dir_var=${wd}/variant
dir_pop=${dir_var}/aory.merge
refgenome=${dir_bam}/AoH_genomic.fa
sample=aory.merge

echo "Aligning ${sample}"
mafft --thread 38 --auto ${dir_pop}/${sample}.fa* > ${dir_pop}/${sample}.align.fa

#build tree from variable data
#GTR General time reversible model with unequal rates and unequal base freq.
#ASC ascertainment bias correction (+ASC) model. +ASC will correct the likelihood conditioned on variable sites. Without +ASC, the branch lengths might be overestimated.
#bnni Optimize UFBoot trees by NNI on bootstrap alignment
iqtree -s ${dir_pop}/${sample}.align.fa --seqtype DNA --runs 1 -T AUTO -B 1000 -bnni --alrt 1000 -m GTR+ASC
