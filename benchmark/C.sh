#!/bin/bash

for i in 260 270 280 290 300 325 350 375 400 500 600 700; do
	echo $i >> nohup.out
	nohup ./caller /data/victor/Syllable-Query/save16.txt load /data6/victor/query16.vcf resultsC${i}.txt sites $i
done
