head -n2 pairlist | while read i a b;
do
Rscript /mnt/content_93/home/yichun_hml/tools/Rscripts/mummerCoordsDotPlotly-grey.R -i $i.coords -o $i.500 -s -t -m 500 -q 500000 -l -p 10 -a $a -b $b
Rscript /mnt/content_93/home/yichun_hml/tools/Rscripts/mummerCoordsDotPlotly-grey.R -i $i.coords -o $i.1000 -s -t -m 1000 -q 500000 -l -p 10 -a $a -b $b
Rscript /mnt/content_93/home/yichun_hml/tools/Rscripts/mummerCoordsDotPlotly-grey.R -i $i.coords -o $i.5000 -s -t -m 5000 -q 500000 -l -p 10 -a $a -b $b
Rscript /mnt/content_93/home/yichun_hml/tools/Rscripts/mummerCoordsDotPlotly-grey.R -i $i.coords -o $i.10000 -s -t -m 10000 -q 500000 -l -p 10 -a $a -b $b
Rscript /mnt/content_93/home/yichun_hml/tools/Rscripts/mummerCoordsDotPlotly-grey.R -i $i.coords -o $i.50000 -s -t -m 50000 -q 500000 -l -p 10 -a $a -b $b
done

##pairlist
AoH-AoRIB40	AoH	AoRIB40
AoL-AoRIB40	AoL	AoRIB40
AoL-AoH	AoL	AoH
AoH-AoL	AoH AoL
AoRIB40-AoH	AoRIB40	AoH
AoRIB40-AoL	AoRIB40	AoL
AoH100x-AoL100x	AoH100x	AoL100x
AoRIB40-AoH100x	AoRIB40	AoH100x
AoRIB40-AoL100x	AoRIB40	AoL100x
AoH10k-AoL10k	AoH10k	AoL10k
AoRIB40-AoH10k	AoRIB40	AoH10k
AoRIB40-AoL10k	AoRIB40	AoL10k
