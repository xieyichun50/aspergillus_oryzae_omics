conda activate cgat-apps

cat namelist3 | while read i;
do 
perl /mnt/content_176/yichun/scripts/comparative/1ortholog_prep/gff2gtf.pl $i*.gff3 | grep -w CDS | sed 's/gene_id \".*\"\; transcript_id/gene_id/' > $i.CDS
cgat gtf2gtf -I $i.CDS  --method=filter --filter-method=longest-gene > $i.longest-gene.gtf
grep -w "gene_id" $i.longest-gene.gtf |sed 's/.*gene_id "//g;s/".*//' |sort -u > $i.longest-gene.gtf.id
seqtk subseq $i_proteins.fa $i.longest-gene.gtf.id > $i.longest-gene.proteins.fa
wc -l $i.longest-gene.gtf.id
wc -l $i.longest-gene.proteins.fa
done

conda deactivate

conda activate orthofinder

orthofinder -f /mnt/content_176/yichun/fungi/aspergillus/comparative -t 76 -a 76 -M msa

conda deactivate

rsync -avzP yichun_hml@137.189.43.176:/mnt/content_176/yichun/fungi/aspergillus/comparative --exclude Gene_Trees/ --exclude MultipleSequenceAlignments/ --exclude Orthogroup_Sequences/ --exclude Resolved_Gene_Trees/ --exclude Single_Copy_Orthologue_Sequences/ --exclude WorkingDirectory/ .
