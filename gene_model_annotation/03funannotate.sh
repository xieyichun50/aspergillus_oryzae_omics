cat genome_list | while read spc;
do
echo "funannotate iprscan -i ${spc} -m local --cpus 72"
funannotate iprscan -i ${spc} -m local --cpus 72 
#echo "interproscan.sh -i ${spc}/annotate_misc/genome.proteins.fasta -f XML -goterms -pa -cpu 72"
#interproscan.sh -i ${spc}/annotate_misc/genome.proteins.fasta -f XML -goterms -pa -cpu 72
done

cat genome_list | while read spc;
do
echo "funannotate annotate -i ${spc} --busco_db ascomycota --cpus 72"
funannotate annotate -i ${spc} --busco_db ascomycota --cpus 72
done
