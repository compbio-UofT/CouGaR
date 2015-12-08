#!/bin/bash

g=/filer/misko/mini_chr/git/minichr/
p=/filer/misko/picard/picard-tools-1.56/
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/
c=/data/misko/2013.04.12/cs2-4.3/cs2.exei

if [ $# -ne 2 ]; then
	echo $0 folder id
	exit
fi

f=$1
id=$2

echo $id `python $g/getcov/plot.py $f/normal_cov.gz.log $f/tumor_cov.gz.log $f/gc.png` > $f/plot.log
echo $id `pypy $g/walker/results_summary.py $f/Qg_ip*292*_m3` > $f/summary.log
