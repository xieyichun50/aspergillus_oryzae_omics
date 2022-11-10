show-snps -CTlr AoH-AoL.delta > AoH-AoL.delta.snp.txt

/mnt/content_93/home/yichun_hml/tools/bin/mummer2vcf.py -s AoH-AoL.delta.snp.txt -g /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/update/AoH/AoH.fasta --no-Ns --type SNP --input-header > AoH-AoL.delta.snp.vcf
java -jar /mnt/content_93/home/yichun_hml/tools/snpEff/snpEff.jar AoH AoH-AoL.delta.snp.vcf > AoH-AoL.delta.snp.anno.vcf

/mnt/content_93/home/yichun_hml/tools/bin/mummer2vcf.py -s AoH-AoL.delta.snp.txt -g /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/update/AoH/AoH.fasta --no-Ns --type INDEL --input-header > AoH-AoL.delta.indel.vcf
java -jar /mnt/content_93/home/yichun_hml/tools/snpEff/snpEff.jar AoH AoH-AoL.delta.indel.vcf > AoH-AoL.delta.indel.anno.vcf

show-snps -CTlr AoL-AoH.delta > AoL-AoH.delta.snp.txt

/mnt/content_93/home/yichun_hml/tools/bin/mummer2vcf.py -s AoL-AoH.delta.snp.txt -g /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/update/AoL/AoL.fasta --no-Ns --type SNP --input-header > AoL-AoH.delta.snp.vcf
java -jar /mnt/content_93/home/yichun_hml/tools/snpEff/snpEff.jar AoL AoL-AoH.delta.snp.vcf > AoL-AoH.delta.snp.anno.vcf

/mnt/content_93/home/yichun_hml/tools/bin/mummer2vcf.py -s AoL-AoH.delta.snp.txt -g /mnt/content_93/home/yichun_hml/aspergillus/genome_assembly/update/AoL/AoL.fasta --no-Ns --type INDEL --input-header > AoL-AoH.delta.indel.vcf
java -jar /mnt/content_93/home/yichun_hml/tools/snpEff/snpEff.jar AoL AoL-AoH.delta.indel.vcf > AoL-AoH.delta.indel.anno.vcf


show-snps -CTlr AoH-AoL.delta > AoH-AoL.delta.snp.txt

/mnt/content_93/home/yichun_hml/tools/bin/mummer2vcf.py -s AoH-AoL.delta.snp.txt -g /mnt/content_93/home/yichun_hml/aspergillus/rawdata/TGS/02.Assembly/AoH/01.assembly/AoH.assembly.seq --no-Ns --type SNP --input-header > AoH-AoL.delta.snp.vcf
java -jar /mnt/content_176/yichun/tools/snpEff/snpEff.jar Aspergillus_oryzae_AoH AoH-AoL.delta.snp.vcf > AoH-AoL.delta.snp.anno.vcf

/mnt/content_93/home/yichun_hml/tools/bin/mummer2vcf.py -s AoH-AoL.delta.snp.txt -g /mnt/content_93/home/yichun_hml/aspergillus/rawdata/TGS/02.Assembly/AoH/01.assembly/AoH.assembly.seq --no-Ns --type INDEL --input-header > AoH-AoL.delta.indel.vcf
java -jar /mnt/content_176/yichun/tools/snpEff/snpEff.jar Aspergillus_oryzae_AoH AoH-AoL.delta.indel.vcf > AoH-AoL.delta.indel.anno.vcf

show-snps -CTlr AoL-AoH.delta > AoL-AoH.delta.snp.txt
/mnt/content_93/home/yichun_hml/tools/bin/mummer2vcf.py -s AoL-AoH.delta.snp.txt -g /mnt/content_93/home/yichun_hml/aspergillus/rawdata/TGS/02.Assembly/AoL/01.assembly/AoL.assembly.seq --no-Ns --type SNP --input-header > AoL-AoH.delta.snp.vcf
java -jar /mnt/content_176/yichun/tools/snpEff/snpEff.jar Aspergillus_oryzae_AoL AoL-AoH.delta.snp.vcf > AoL-AoH.delta.snp.anno.vcf

/mnt/content_93/home/yichun_hml/tools/bin/mummer2vcf.py -s AoL-AoH.delta.snp.txt -g /mnt/content_93/home/yichun_hml/aspergillus/rawdata/TGS/02.Assembly/AoL/01.assembly/AoL.assembly.seq --no-Ns --type INDEL --input-header > AoL-AoH.delta.indel.vcf
java -jar /mnt/content_176/yichun/tools/snpEff/snpEff.jar Aspergillus_oryzae_AoL AoL-AoH.delta.indel.vcf > AoL-AoH.delta.indel.anno.vcf
