#!/bin/bash

Ls=(936 1337 1571 1791 1926 2078)

# for i in 0 1 2 3 4 5; do
for i in 1 2; do
	echo $i >> nohup_fullN${i}.out
	nohup ./d-PBWT3 panel-${i}.vcf query-${i}.vcf ${Ls[$i]} N_fullres${i}.txt >> nohup_fullN${i}.out 2>&1
done
