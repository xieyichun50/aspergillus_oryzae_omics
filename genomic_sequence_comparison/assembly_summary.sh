## assembly-stats
ls *.fasta | while read i;
do
assembly-stats -t $i
done

## busco
cat /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/filelist | while read i;
do 
busco -i /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/$i.fasta -o busco.ascomycota_odb10.$i -m genome -l ascomycota_odb10 --download_path /mnt/content_93/home/yichun_hml/tools/busco_db --augustus_species aspergillus_oryzae -c 12 --offline
busco -i /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/$i.fasta -o busco.eukaryota_odb10.$i -m genome -l eukaryota_odb10 --download_path /mnt/content_93/home/yichun_hml/tools/busco_db --augustus_species aspergillus_oryzae -c 12 --offline
busco -i /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/$i.fasta -o busco.eurotiales_odb10.$i -m genome -l eurotiales_odb10 --download_path /mnt/content_93/home/yichun_hml/tools/busco_db --augustus_species aspergillus_oryzae -c 12 --offline
done

