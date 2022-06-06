#!/bin/bash

for i in 0 1 2 3 4 5 6 7; do
	echo $i >> nohup_128M${i}.out
	nohup ./caller128 panelM-${i}.vcf queryM-${i}.vcf M_128res${i}.txt 2000 >> nohup_128M${i}.out 2>&1
done
