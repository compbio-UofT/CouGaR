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

cov_mapq=20

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
	dname=`echo $line  | sed 's/.*download[/]\(.*\)/\1/g'`
	#echo "dname is $dname"
	$c "analysis_id=${line}" -o get.xml
	$gt -k 5 --peer-timeout 30 -v -c $key -d get.xml
	bamfile=$dname/*.bam
	ls -l $bamfile
	#check the MD5
	echo $bamfile
	if [ ! -e $bamfile ] ; then 
		echo "DID NOT GET BAM FILE!"
		rm normal.bam
		rm tumor.bam
		rm ${dwndir}/$tid
		touch fail
		exit 100 #SGE ERROR
	fi
	md5=`md5sum $bamfile | awk '{print $1}'`
	echo MD5 is $md5
	echo "XXXXX"
	grep -i md5 get.xml
	echo "XXXXX"
	grep -i md5 downloadable.xml
	echo "XXXXX"
	grep -i $md5 get.xml > /dev/null
	if [ $? -eq 0 ]; then 
		echo "md5 is ok"
	else 
		echo "md5 is not ok"
		rm normal.bam
		rm tumor.bam
		rm ${dwndir}/$tid
		touch fail
		exit 100 #SGE ERROR
	fi
	rm get.xml
	echo $dname, $bamfile
	if [ ! -e $bamfile ] ; then 
		echo DID NOT DOWNLOAD BAM FILE! ERROR
		rm normal.bam
		rm tumor.bam
		rm ${dwndir}/$tid
		touch fail
		exit 100 #SGE ERROR
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
		$s mpileup -q ${cov_mapq} tumor.bam | $g/getcov/get_cov tumor_cov 
		#gzip tumor_cov
		sh $g/make_clusters.sh tumor.bam 0 
		sh $g/make_clusters.sh tumor.bam 15
		rm tumor.bam
	else
		echo $bamfile normal
		mv $bamfile normal.bam
		$s index normal.bam
		$s mpileup -q ${cov_mapq} normal.bam | $g/getcov/get_cov normal_cov 
		#gzip normal_cov
		sh $g/make_clusters.sh normal.bam 0 
		sh $g/make_clusters.sh normal.bam 15 
		rm normal.bam
	fi
done

#remove the download lock
rm ${dwndir}/$tid

