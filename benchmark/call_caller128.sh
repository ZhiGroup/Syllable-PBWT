#!/bin/bash

# for i in 260 270 280 290 300 325 350 375 400 500 600 700; do
for i in 270 290; do
	echo $i >> nohup_uk16_128.out
	nohup ./caller128 uk16_panel.vcf uk16_query.vcf uk16_128res${i}.txt $i >> nohup_uk16_128.out 2>&1
done
