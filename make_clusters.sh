#!/bin/bash

g=/dupa-filer/misko/cougar/CouGaR/
p=/filer/misko/picard/picard-tools-1.56/
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/
c=/data/misko/2013.04.12/cs2-4.3/cs2.exe

#using hg19
#ref=/filer/hg19/hg19.fa
#or hg18
#ref=/filer/hg18/hg18.fa

if [ $# -ne 2 ]; then 
	echo $0 bam clustermapq
	exit
fi

cluster_mapq=$6


bam=$1
if [ ! -f $bam ]; then
	echo Cannot find bam file!
	#exit
fi
if [ ! -f $bam.bai ]; then
	echo Bam file not indexed!
	#exit
fi
clustermapq=$2


cluster=$g/clustering/cluster

echo `date` $0 $@ " making clusters clustermapq:" ${clustermapq} >> ${bam}_command_line

function mcluster {
	msfile=`echo $bam | sed 's/.bam/_mean_and_std.txt/g'`
	if [ ! -f $msfile -o `cat $msfile | awk 'END {print NR}'` -lt 1 ]; then
		echo failed to find mean and std file!
		$j -jar $p/CollectInsertSizeMetrics.jar VALIDATION_STRINGENCY=SILENT H=${bam}_histo I=${bam} O=${bam}_stats AS=true 
		cat ${bam}_stats | awk '{print $5,$6}' | grep -A 1 DEVI | tail -n 1 > ${msfile}
	fi
	if [ ! -f $msfile ] ; then
		echo failed to make insert size metrics
		exit
	fi
	local mean=`cat $msfile | awk '{print int($1)}'`
	local stddev=`cat $msfile | awk '{print int($2)}'`
	
	rfn=`echo $bam | sed 's/.bam//g'`

	outd=${rfn}_clusters
	mkdir -p $outd
	outa=${rfn}_clusters/q${clustermapq}.txt.gz
	if [ -e $outa ]; then
	sz=`/usr/bin/du -s $outa | awk '{print $1}'`
	else
	sz=0
	fi
	if [ ! -e $outa -o $sz -lt 30 ]; then 
		rm $outa
		echo Generating base clusters
		$s view -q ${clustermapq} ${bam} | python $g/filter_bwa_mq.py ${clustermapq} | $cluster $mean $stddev | gzip > ${outa}_wchrM
		zgrep -v chrM ${outa}_wchrM | gzip > $outa
	else
		echo Skipping gen
	fi
	covs="10 5 3 2 0"
	for cov in $covs; do
		#move over the clusters, filter for same strand (type 0 or type 1) with size less then 2000
		out=${rfn}_clusters/q${clustermapq}_cov${cov}.txt.gz
		zcat $outa | awk -v c=$cov '{if ($4>c) {print $0}}' | sed 's/\([0-9]\)\t\([0-9]*\)[:]\([0-9]*\)\t\([0-9]*\)[:]\([0-9]*\)\t\([0-9]*\)\t\(.*\)\t\([0-9][0-9]*[.]*[0-9]*\)\t\([0-9][0-9]*[.]*[0-9]*\)/\2\t\3\t\4\t\5\t\1\t\6\t\8\t\9/g' | awk '{if ($1==23) {$1="X"}; if ($1==24) {$1="Y"}; if ($3==23) {$3="X"}; if ($3==24) {$3="Y"}; a="+"; b="+"; if ($5==1) {a="-"; b="-"}; if ($5==2) {a="+"; b="-"}; if ($5==3) {a="-"; b="+"}; print "chr"$1"\t"$2"\t"a"\tchr"$3"\t"$4"\t"b"\t"$6"\t"$7"\t"$8}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}' | gzip > $out
		zcat $out | awk '{if ($8<200 && $9<200) {print $0}}' | gzip > ${out}_200.gz	
	done
}



##generate the rough clusters
echo current dir is `pwd`
mcluster 
