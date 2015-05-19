#!/bin/bash


c=/usr/bin/cgquery
gt=/usr/bin/gtdownload
key=/dupa-filer/misko/tcga/cghub.key 
s=/home/misko/apps/bin/samtools
g=/filer/misko/mini_chr/git/minichr

if [ $# -ne 2 ]; then
	echo $0 TCGA_id wd
	exit
fi

tid=$1
wd=$2


mkdir -p $wd 
pushd $wd


#grab a download space
dwndir=/dupa-filer/misko/tcga_data/downloads/



if [ -e normal.bam -a -e tumor.bam ]; then
	echo found both!
fi

cov_mapq=20
$s mpileup -q ${cov_mapq} tumor.bam | $g/getcov/get_cov tumor_cov &
$s mpileup -q ${cov_mapq} normal.bam | $g/getcov/get_cov normal_cov &
wait

echo Generating clusters...
sh $g/make_clusters.sh tumor.bam 0 
sh $g/make_clusters.sh tumor.bam 15
sh $g/make_clusters.sh normal.bam 0 
sh $g/make_clusters.sh normal.bam 15 
	

#rm normal.bam
#rm tumor.bam
