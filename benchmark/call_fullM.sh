#!/bin/bash

for i in 0 1 2 3; do
	echo $i >> nohup_fullM${i}.out
	nohup ./d-PBWT3 panelM-${i}.vcf queryM-${i}.vcf 2000 M_fullres${i}.txt >> nohup_fullM${i}.out 2>&1
done
