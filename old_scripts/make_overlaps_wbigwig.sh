#!/bin/bash

g=/filer/misko/mini_chr/git/minichr/
p=/filer/misko/picard/picard-tools-1.56/
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/
c=/data/misko/2013.04.12/cs2-4.3/cs2.exe

#using hg19
#ref=/filer/hg19/hg19.fa
#or hg18
#ref=/filer/hg18/hg18.fa

if [ $# -ne 1 ]; then 
	echo $0 folder
	exit
fi
wd=$1

function moverlap {
	pushd $wd
	python $g/scripts/overlap.py 3000 ${wd}/normal_clusters/q0_cov0.txt.gz ${wd}/tumor_clusters/q0_cov5.txt.gz_200.gz  |  awk '{if (NF>4) {print $0}}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}'  > $wd/nsubtract_centrosubtract_1000bp
	cat $wd/nsubtract_centrosubtract_1000bp | sed 's/\(chr[^:]*\):\([0-9]*\)\([+-]\)\s\(chr[^:]*\):\([0-9]*\)\([+-]\)/\1\t\2\t\3\t\4\t\5\t\6/g' | awk '{OFS="\t"; type=0; if ($3=="+") {if ($6=="+") {type=0} else {type=2} } else { if ($6=="+") {type=3} else {type=1} }; print $1,$2,$5,type,$7,0,0,0.0,0,$4,"EDGE" }' > $wd/nsubtract_centrosubtract_1000bp_links
	$g/hmm $wd/nsubtract_centrosubtract_1000bp_links ${wd}/tumor_cov.gz ${wd}/normal_cov.gz 0 > ${wd}/hmm

	#compute the edge mapqs
        while read line; do 
                #s=`echo $line | awk '{if ($3=="+") {print $1"\t"$2-100"\t"$2+20} else {print $1"\t"$2-20"\t"$2+100}}' | sed 's/ /\t/g'`
                #e=`echo $line | awk '{if ($6=="+") {print $4"\t"$5-20"\t"$5+100} else {print $4"\t"$5-100"\t"$5+20}}' | sed 's/ /\t/g'`
		echo `echo $line | awk '{s=$2-30; if (s<0) {s=1}; print $1,s,$2+30,$1":"$2}'` 
		echo `echo $line | awk '{s=$5-30; if (s<0) {s=1}; print $4,s,$5+30,$4":"$5}'` 
		#echo $e `echo $line | awk '{print $4":"$5$6}'`
                #echo `echo $line | awk '{print $1":"$2}'` `$s view tumor.bam $s | awk '{s+=$5} END {print s/NR}'`
                #echo `echo $line | awk '{print $4":"$5}'` `$s view tumor.bam $e | awk '{s+=$5} END {print s/NR}'` 
        done < nsubtract_centrosubtract_1000bp  | sort | uniq  > nsubtract_centrosubtract_1000bp.bed
	ref=`python $g/get_ref.py downloadable.xml | awk '{print $NF}' | head -n 1`
	if [ -z "$ref" ] ;  then
		echo "FAILED!!! COULD NOT FIND REF"
		exit
	fi
	$g/mapability/bigWigAverageOverBed $g/mapability/$ref/* nsubtract_centrosubtract_1000bp.bed nsubtract_centrosubtract_1000bp.bed.out
	cat nsubtract_centrosubtract_1000bp.bed.out | awk '{print $1,$NF}' > nsubtract_centrosubtract_1000bp_mqs 


	$g/getcov/print_cov ${wd}/tumor_cov.gz > ${wd}/tumor_cov.gz.log
	$g/getcov/print_cov ${wd}/normal_cov.gz > ${wd}/normal_cov.gz.log
	# | sort | uniq > nsubtract_centrosubtract_1000bp_mqs
	popd
}

#run overlaps
moverlap




