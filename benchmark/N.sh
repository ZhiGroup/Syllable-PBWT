#!/bin/bash

Ls=(0 936 1337 1571 1791 1926 2078)
for i in {1..6}; do
	nohup ./caller /data6/victor/panelN${i}.vcf saveN${i}.txt /data6/victor/queryN${i}.vcf resultsN${i}.txt sites ${Ls[$i]}
done
