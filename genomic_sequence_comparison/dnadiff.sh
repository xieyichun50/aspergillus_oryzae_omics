##no filter
dnadiff -p AoH-AoL /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.fasta
show-coords -c AoH-AoL.delta > AoH-AoL.coords

dnadiff -p AoRIB40-AoH /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/AoRIB40.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.fasta
show-coords -c AoRIB40-AoH.delta > AoRIB40-AoH.coords

dnadiff -p AoRIB40-AoL /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/AoRIB40.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.fasta
show-coords -c AoRIB40-AoL.delta > AoRIB40-AoL.coords

dnadiff -p AoL-AoH /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.fasta
show-coords -c AoL-AoH.delta > AoL-AoH.coords

dnadiff -p AoH-AoRIB40 /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/AoRIB40.fasta
show-coords -c AoH-AoRIB40.delta > AoH-AoRIB40.coords

dnadiff -p AoL-AoRIB40 /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/AoRIB40.fasta
show-coords -c AoL-AoRIB40.delta > AoL-AoRIB40.coords

##Longest 100x
dnadiff -p AoH100x-AoL100x /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.100x.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.100x.fasta
show-coords -c AoH100x-AoL100x.delta > AoH100x-AoL100x.coords

dnadiff -p AoRIB40-AoH100x /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/AoRIB40.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.100x.fasta
show-coords -c AoRIB40-AoH100x.delta > AoRIB40-AoH100x.coords

dnadiff -p AoRIB40-AoL100x /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/AoRIB40.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.100x.fasta
show-coords -c AoRIB40-AoL100x.delta > AoRIB40-AoL100x.coords

##Over 10k reads
dnadiff -p AoH10k-AoL10k /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.10k.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.10k.fasta
show-coords -c AoH10k-AoL10k.delta > AoH10k-AoL10k.coords

dnadiff -p AoRIB40-AoH10k /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/AoRIB40.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoH.10k.fasta
show-coords -c AoRIB40-AoH10k.delta > AoRIB40-AoH10k.coords

dnadiff -p AoRIB40-AoL10k /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/AoRIB40.fasta /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/TGS.AoL.10k.fasta 
show-coords -c AoRIB40-AoL10k.delta > AoRIB40-AoL10k.coords
