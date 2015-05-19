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

if [ -d $wd -a `ls $wd | wc -l | awk '{print $1}'` -ne 1 ]; then 
	echo Directory exists! $wd
	exit
fi

mkdir -p $wd 
pushd $wd


$c "filename=*${tid}*&library_strategy=WGS" -o downloadable.xml
files=`cat downloadable.xml | grep downloadable_file_count | sed 's/.*>\([0-9]*\)<.*/\1/g'`
echo There are $files files to be downloaded 
# if [ $files -ne 4 ] then 
#	echo "File downloaded aborted wrong count"
#	exit
# fi


found_normal=0
found_tumor=0

#grab a download space
dwndir=/dupa-filer/misko/tcga_data/downloads/

candownload=0

echo Grabbing download lock...
#grab a download lock
while [ $candownload -eq 0 ]; do
	max_downloads=`cat ${dwndir}/max_downloads`
	sleep $[ ( $RANDOM % 60 )  + 1 ]
	n=`ls ${dwndir} | wc -l`
	echo waiting for download lock $n
	if [ $n -lt ${max_downloads} ] ; then
		touch ${dwndir}/$tid
		candownload=1
	fi
done
echo Have download lock...


sp=`/bin/df -m /dupa-filer/  | grep dupa | awk '{print $3}'`
while [ $sp -lt 1200000 ] ; do
	sleep 60
	echo waiting for space
	sp=`/bin/df -m /dupa-filer/  | grep dupa | awk '{print $3}'`
done

python /filer/misko/mini_chr/git/minichr/get_analysis_ids.py downloadable.xml | while read line; do 
	$c "analysis_id=${line}" -o get.xml
	$gt -v -c $key -d get.xml
	rm get.xml
	dname=`echo $line  | sed 's/.*download[/]\(.*\)/\1/g'`
	bamfile=$dname/*.bam
	echo $dname, $bamfile
	if [ ! -e $bamfile ] ; then 
		echo DID NOT DOWNLOAD BAM FILE! ERROR
		rm normal.bam
		rm tumor.bam
		rm ${dwndir}/$tid
		exit
	fi
	info=`echo $bamfile | sed 's/.*TCGA-\([A-Z0-9][A-Z0-9]\)-\([A-Z0-9][A-Z0-9][A-Z0-9][A-Z0-9]\)-\([A-Z0-9][A-Z0-9]\)\([A-Z]\)-.*/\1 \2 \3 \4/g'`
	tss=`echo $info | awk '{print $1}'`
	participant=`echo $info | awk '{print $2}'`
	sample=`echo $info | awk '{print $3}'`
	vial=`echo $info | awk '{print $4}'`
	echo $tss X $participant X $sample X $vial X $info

	if [ $sample -lt 10 ] ; then
		echo $bamfile tumor
		mv $bamfile tumor.bam
		$s index tumor.bam
		found_tumor=1
	else
		echo $bamfile normal
		mv $bamfile normal.bam
		$s index normal.bam
		found_normal=1
	fi
done

#remove the download lock
rm ${dwndir}/$tid

if [ -e normal.bam -a -e tumor.bam ]; then
	echo found both!
fi
cov_mapq=20
$s mpileup -q ${cov_mapq} tumor.bam | $g/getcov/get_cov tumor_cov &
$s mpileup -q ${cov_mapq} normal.bam | $g/getcov/get_cov normal_cov &
echo Generating coverages...
wait

echo Generating clusters...
sh $g/make_clusters.sh tumor.bam 0 
sh $g/make_clusters.sh tumor.bam 15
sh $g/make_clusters.sh normal.bam 0 
sh $g/make_clusters.sh normal.bam 15 
	

rm normal.bam
rm tumor.bam
