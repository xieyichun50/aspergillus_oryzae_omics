## snpEff installation note
```
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

In `snpEff.config`
```
#Aspergillus_oryzae_AoH
Aspergillus_oryzae_AoH.genome:Aspergillus_oryzae_AoH

#Aspergillus_oryzae_AoL
Aspergillus_oryzae_AoL.genome:Aspergillus_oryzae_AoL

```
```
mkdir data
cd data
mkdir genomes
mkdir Aspergillus_oryzae_AoH
mkdir Aspergillus_oryzae_AoL
rsync /mnt/content_176/yichun/fungi/aspergillus/genome/AoH_genomic.fa genomes/Aspergillus_oryzae_AoH.fa
rsync /mnt/content_176/yichun/fungi/aspergillus/genome/AoL_genomic.fa genomes/Aspergillus_oryzae_AoL.fa
rsync /mnt/content_176/yichun/fungi/aspergillus/genome/AoH_genomic.gff3 Aspergillus_oryzae_AoH/genes.gff
rsync /mnt/content_176/yichun/fungi/aspergillus/genome/AoH_proteins.fa  Aspergillus_oryzae_AoH/protein.fa
rsync /mnt/content_176/yichun/fungi/aspergillus/genome/AoL_genomic.gff3 Aspergillus_oryzae_AoL/genes.gff
rsync /mnt/content_176/yichun/fungi/aspergillus/genome/AoL_proteins.fa  Aspergillus_oryzae_AoL/protein.fa
```
```
java -jar /mnt/content_176/yichun/tools/snpEff/snpEff.jar build -gff3 -v Aspergillus_oryzae_AoH
java -jar /mnt/content_176/yichun/tools/snpEff/snpEff.jar build -gff3 -v Aspergillus_oryzae_AoL
```

In the `$dir`

```
/mnt/content_176/yichun/tools/snpEff/data/
├── Aspergillus_oryzae_AoH
│   ├── genes.gff
│   ├── protein.fa
│   └── snpEffectPredictor.bin
├── Aspergillus_oryzae_AoL
│   ├── genes.gff
│   ├── protein.fa
│   └── snpEffectPredictor.bin
└── genomes
    ├── Aspergillus_oryzae_AoH.fa
    └── Aspergillus_oryzae_AoL.fa

3 directories, 8 files

```
