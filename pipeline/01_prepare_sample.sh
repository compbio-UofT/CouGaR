#!/bin/bash

#check that COUGARD has been properly set
if [ -z "$COUGARD" ]; then echo "Please set COUGAR directory COUGARD=path"; exit; fi
if [ ! -d "$COUGARD" ] ; then echo "COUGARD enviornment variable is not set properly"; exit; fi

#load the configuation 
. $COUGARD/cougar_conf.sh



if [ $# -ne 5 ]; then
	echo $0 "ID_folder_name tumor_bam normal_bam ref[hg18/hg19] group[anygrouplabel]"
	exit
fi

wd=$1
tumor_bamfilename=$2
normal_bamfilename=$3
ref=$4
group=$5

if [ "$ref" == "hg18" ] ; then 
	echo "is hg18"
elif [ "$ref" == "hg19" ] ; then 
	echo "is hg19"
else
	echo "Reference genome is unrecognized, please use either hg19 or hg18"
	exit
fi

if [ -d $wd -a `ls $wd | wc -l | awk '{print $1}'` -ne 1 ]; then 
	echo Directory exists! $wd
	#exit
fi

mkdir -p $wd 
pushd $wd
	echo $ref > ref
	echo $group > subset
	bamfile=$tumor_bamfilename
	echo $bamfile tumor
	ln -s $bamfile tumor.bam
	$s index tumor.bam
	$s mpileup -q ${cov_mapq} tumor.bam | $g/getcov/get_cov tumor_cov 
	sh $g/make_clusters.sh tumor.bam 0 &
	sh $g/make_clusters.sh tumor.bam 15 &

	bamfile=$normal_bamfilename
	echo $bamfile normal
	ln -s $bamfile normal.bam
	$s index normal.bam
	$s mpileup -q ${cov_mapq} normal.bam | $g/getcov/get_cov normal_cov 
	sh $g/make_clusters.sh normal.bam 0 &
	sh $g/make_clusters.sh normal.bam 15 &
wait
	rm tumor.bam
	rm normal.bam
popd
