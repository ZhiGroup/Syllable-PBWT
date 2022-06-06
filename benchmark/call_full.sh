#!/bin/bash

# for i in 260 270 280 290 300 325 350 375 400 500 600 700; do
for i in 700; do
	echo $i >> nohup_uk16_full.out
	nohup ./d-PBWT3 uk16_panel.vcf uk16_query.vcf $i uk16_fullres${i}.txt >> nohup_uk16_full.out 2>&1
done
