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

cov_mapq=20

if [ "$ref" != "hg18" -a "$ref" != "hg19" ] ; then 
	echo "Reference genome is unrecognized, please use either hg19 or hg18"
	exit
fi

if [ -d "$wd" ]; then 
	echo Directory exists! $wd
	#exit
fi

if [ ! -e ${tumor_bamfilename} ] ; then 
	echo Tumor bam file \"${tumor_bamfilename}\" does not exist
	exit
fi

if [ ! -e ${normal_bamfilename} ] ; then 
	echo normal bam file \"${normal_bamfilename}\" does not exist
	exit
fi

mkdir -p $wd 
pushd $wd
	echo $ref > ref
	echo $group > subset
	bamfile=$tumor_bamfilename
	echo Using $bamfile as tumor BAM ...
	ln -s $bamfile tumor.bam
	echo Indexing tumor BAM
	$s index tumor.bam
	echo Creating coverage file for tumor BAM
	$s mpileup -q ${cov_mapq} tumor.bam 2>/dev/null | $g/getcov/get_orig_cov tumor_orig_cov 
	echo Background cluster generation for tumor BAM
	sh $g/clustering/make_clusters.sh tumor.bam 0 &

	bamfile=$normal_bamfilename
	echo Using $bamfile as normal BAM ...
	ln -s $bamfile normal.bam
	echo Indexing normal BAM
	$s index normal.bam
	echo Creating coverage file for normal BAM
	$s mpileup -q ${cov_mapq} normal.bam 2>/dev/null | $g/getcov/get_cov normal_cov 
	echo Background cluster generation for normal BAM
	sh $g/clustering/make_clusters.sh normal.bam 0 &
echo Waiting for cluster generation to finish ... 
wait
echo All done
	#rm tumor.bam
	#rm normal.bam
popd
