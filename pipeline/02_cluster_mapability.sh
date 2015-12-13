#!/bin/bash

#check that COUGARD has been properly set
if [ -z "$COUGARD" ]; then echo "Please set COUGAR directory COUGARD=path"; exit; fi
if [ ! -d "$COUGARD" ] ; then echo "COUGARD enviornment variable is not set properly"; exit; fi

#load the configuation 
. $COUGARD/cougar_conf.sh


if [ $# -ne 1 ]; then 
	echo $0 folder
	exit
fi
wd=$1

if [ ! -d $wd ] ; then 
	echo Directoty $wd does not exist!
	exit
fi

function moverlap {
	pushd $wd

	ref=`cat ref | head -n 1`
	if [ -z "$ref" ] ;  then
		echo "FAILED!!! COULD NOT FIND REF"
		exit
	fi

	

	#find clusters that are not in the normal sample but only in the tumor, approximately
	python $g/scripts/overlap.py 3000 normal_clusters/q0_cov0.txt.gz tumor_clusters/q0_cov5.txt.gz_200.gz  |  awk '{if (NF>4) {print $0}}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}'  > nsubtract_centrosubtract_1000bp

        #convert the above results to a different format
	cat nsubtract_centrosubtract_1000bp | sed 's/\(chr[^:]*\):\([0-9]*\)\([+-]\)\s\(chr[^:]*\):\([0-9]*\)\([+-]\)/\1\t\2\t\3\t\4\t\5\t\6/g' | awk '{OFS="\t"; type=0; if ($3=="+") {if ($6=="+") {type=0} else {type=2} } else { if ($6=="+") {type=3} else {type=1} }; print $1,$2,$5,type,$7,0,0,0.0,0,$4,"EDGE" }' > nsubtract_centrosubtract_1000bp_links

	#run the hmm using the tumor and normal coverage by also applying priors for copy number transistion on cluster break point locations
	$g/hmm/hmm_chrm nsubtract_centrosubtract_1000bp_links tumor_cov.gz normal_cov.gz 0 $ref > hmm

	#compute mapqs around the breakpoints of each cluster
        while read line; do 
		echo `echo $line | awk '{s=$2-30; if (s<0) {s=1}; print $1,s,$2+30,$1":"$2}'` 
		echo `echo $line | awk '{s=$5-30; if (s<0) {s=1}; print $4,s,$5+30,$4":"$5}'` 
        done < nsubtract_centrosubtract_1000bp  | sort | uniq  > nsubtract_centrosubtract_1000bp.bed
	#ref=`python $g/get_ref.py downloadable.xml | awk '{print $NF}' | head -n 1`
	$g/mapability/bigWigAverageOverBed $g/mapability/$ref/* nsubtract_centrosubtract_1000bp.bed nsubtract_centrosubtract_1000bp.bed.out
	cat nsubtract_centrosubtract_1000bp.bed.out | awk '{print $1,$NF}' > nsubtract_centrosubtract_1000bp_mqs 

	#create a log of coverage for both samples that is human readable
	$g/getcov/print_cov tumor_cov.gz $hgref/${ref}.fa.gz > tumor_cov.gz.log
	$g/getcov/print_cov normal_cov.gz $hgref/${ref}.fa.gz > normal_cov.gz.log
	popd
}

#run overlaps
moverlap
