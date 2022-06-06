#!/bin/bash

Ls=(936 1337 1571 1791 1926 2078)

for i in 0 1 2 3 4 5; do
	echo $i >> nohup_128N${i}.out
	nohup ./caller128 panel-${i}.vcf query-${i}.vcf N_128res${i}.txt ${Ls[$i]} >> nohup_128N${i}.out 2>&1
done
