#!/bin/bash

for i in 260 270 280 290 300 325 350 375 400 500 600 700; do
	echo $i >> nohup_uk16_64.out
	nohup ./caller64 uk16_panel.vcf uk16_query.vcf uk16_64res${i}.txt $i >> nohup_uk16_64.out 2>&1
done
