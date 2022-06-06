#!/bin/bash

for i in 0 1 2 3 4 5 6 7; do
	echo $i >> nohup_64M${i}.out
	nohup ./caller64 panelM-${i}.vcf queryM-${i}.vcf M_64res${i}.txt 2000 >> nohup_64M${i}.out 2>&1
done
