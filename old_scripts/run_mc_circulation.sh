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

if [ $# -ne 6 ]; then 
	echo $0 normal_bam tumor_bam work_directory ref centro clustermapq
	exit
fi

cr=$5
cluster_mapq=$6

ref=$4
if [ ! -f $ref ]; then
	echo Cannot find ref sequence
	exit
fi 
echo using ref $ref
sleep 3

normalbam=$1
if [ ! -f $normalbam ]; then
	echo Cannot find normal bam file!
	#exit
fi
if [ ! -f $normalbam.bai ]; then
	echo Bam file not indexed!
	#exit
fi
tumorbam=$2
if [ ! -f $tumorbam ]; then
	echo Cannot find normal bam file!
	#exit
fi
if [ ! -f $tumorbam.bai ]; then
	echo Bam file not indexed!
	#exit
fi
wd=$3
if [ -d ${wd} ] ; then
	echo Temp directory exists!
	#exit
fi

mkdir -p $wd

#link in the files
if [ ! -e $wd/normal.bam ] ; then 
	rm $wd/normal.bam
	if [ ! -e $normalbam ]; then
		echo error with normal bam
		exit 1
	fi 
	ln -s $normalbam $wd/normal.bam
fi
if [ ! -e $wd/normal.bam.bai ] ; then 
	rm $wd/normal.bam.bai
	if [ ! -e ${normalbam}.bai ]; then
		echo error with normal bai
		$s index $wd/normal.bam
		if [ ! -e $wd/normal.bam.bai ] ; then 
			echo failed to index exiting..
			exit 1
		fi		
	else 
		ln -s ${normalbam}.bai $wd/normal.bam.bai
	fi
fi
if [ ! -e $wd/tumor.bam ] ; then 
	rm $wd/tumor.bam 
	if [ ! -e $tumorbam ]; then
		echo error with tumor bam
		exit 1
	fi 
	ln -s $tumorbam $wd/tumor.bam
fi
if [ ! -e $wd/tumor.bam.bai ] ; then 
	rm $wd/tumor.bam.bai
	if [ ! -e ${tumorbam}.bai ]; then
		echo error with tumor bai
		$s index $wd/tumor.bam
		if [ ! -e $wd/tumor.bam.bai ] ; then 
			echo failed to index exiting..
			exit 1
		fi		
	else 
		ln -s ${tumorbam}.bai $wd/tumor.bam.bai
	fi
fi

tys="normal tumor"
cluster=$g/clustering/cluster_x
cov_mapq=20

echo `date` $0 $@ ",clustermapq:" ${cluster_mapq} ",cov_mapq:" ${cov_mapq} >> $wd/command_line

function covNcluster {
	local ty=$1
	#BAM file coverage
	if [ ! -e ${wd}/${ty}_cov ]; then
		rm ${wd}/${ty}_cov
		echo Generating coverage file
		$s mpileup -q ${cov_mapq} $wd/${ty}.bam | $g/getcov/get_cov ${wd}/${ty}_cov
	else
		echo Skipping cov file generation
	fi

	if [ ! -e $wd/${ty}_stats ] ; then
		rm ${wd}/${ty}_stats
		echo Generating stats for bam
		$j -jar $p/CollectInsertSizeMetrics.jar I=$wd/${ty}.bam H=$wd/${ty}_histo O=$wd/${ty}_stats VALIDATION_STRINGENCY=LENIENT
	else
		echo Skipping stats generation
	fi
	
	cat $wd/${ty}_stats | grep -A 1 MED | grep -v "\-\-" | grep -v MEDIAN | awk '{m+=$5; s+=$6} END {print m/NR,s/NR}' > $wd/${ty}_mean_and_std.txt
	local mean=`cat ${wd}/${ty}_mean_and_std.txt | awk '{print int($1)}'`
	local stddev=`cat ${wd}/${ty}_mean_and_std.txt | awk '{print int($2)}'`
	

	if [ ! -e ${wd}/${ty}_clusters.txt.gz -o ${wd}/${ty}_cov -nt ${wd}/${ty}_clusters.txt.gz ]; then
		rm ${wd}/${ty}_clusters.txt.gz 
		echo Generating base clusters
		$s view -q ${cluster_mapq} ${wd}/${ty}.bam | grep -v chrM  | python $g/filter_bwa_mq.py ${cluster_mapq} | $cluster $mean $stddev | gzip > ${wd}/${ty}_clusters.txt.gz
	else
		echo Skipping base clusters generation
	fi
	#echo SKIPPING MATE MQ FILTER
	#$s view -q ${cluster_mapq} ${wd}/${ty}.bam | grep -v chrM  | $cluster $mean $stddev | gzip > ${wd}/${ty}_clusters.txt.gz
	covs="10 5 3 2 0"
	for cov in $covs; do
		#zcat ${wd}/${ty}_clusters.txt.gz | awk -v c=$cov '{if ($4>c) {print $0}}' | sed 's/\([0-9]\)\t\([0-9]*\)[:]\([0-9]*\)\t\([0-9]*\)[:]\([0-9]*\)\t\([0-9]*\)\t\(.*\)/\2\t\3\t\4\t\5\t\1\t\6/g' | awk '{if ($1==23) {$1="X"}; if ($1==24) {$1="Y"}; if ($3==23) {$3="X"}; if ($3==24) {$3="Y"}; a="+"; b="+"; if ($5==1) {a="-"; b="-"}; if ($5==2) {a="+"; b="-"}; if ($5==3) {a="-"; b="+"}; print "chr"$1"\t"$2"\t"a"\tchr"$3"\t"$4"\t"b"\t"$6}' |  gzip > ${wd}/${ty}_cov${cov}.txt_f.gz
		#move over the clusters, filter for same strand (type 0 or type 1) with size less then 2000
		zcat ${wd}/${ty}_clusters.txt.gz | awk -v c=$cov '{if ($4>c) {print $0}}' | sed 's/\([0-9]\)\t\([0-9]*\)[:]\([0-9]*\)\t\([0-9]*\)[:]\([0-9]*\)\t\([0-9]*\)\t\(.*\)\t\([0-9][0-9]*[.]*[0-9]*\)\t\([0-9][0-9]*[.]*[0-9]*\)/\2\t\3\t\4\t\5\t\1\t\6\t\8\t\9/g' | awk '{if ($1==23) {$1="X"}; if ($1==24) {$1="Y"}; if ($3==23) {$3="X"}; if ($3==24) {$3="Y"}; a="+"; b="+"; if ($5==1) {a="-"; b="-"}; if ($5==2) {a="+"; b="-"}; if ($5==3) {a="-"; b="+"}; print "chr"$1"\t"$2"\t"a"\tchr"$3"\t"$4"\t"b"\t"$6"\t"$7"\t"$8}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}' | gzip > ${wd}/${ty}_cov${cov}.txt_f.gz
	done
}



##generate the rough clusters
for ty in $tys; do
	covNcluster $ty &
done
wait

###python2.6 $g/scripts/overlap.py 1000 ${wd}/normal_cov10.txt_f.gz ${wd}/tumor_cov5.txt_f.gz  |  awk '{if (NF>4) {print $0}}' | awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'  | $b/bin/intersectBed -a - -b $g/centromeres_merged -v | awk '{print $5"\t"$6"\t"$6+1"\t"$7"\t"$1"\t"$2"\t"$4"\t"$8}' | $b/bin/intersectBed -a - -b $g/centromeres_merged -v  | awk '{print $5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$4"\t"$8}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}' > $wd/nsubtract_centrosubtract_1000bp

#echo not removing centromeres

if [ ! -e ${wd}/nsubtract_centrosubtract_1000bp_joined -o ${wd}/${ty}_clusters.txt.gz -nt ${wd}/nsubtract_centrosubtract_1000bp_joined ] ; 
then	
	#take out clusters that have std of 200 or higher in breakpoints (typically centromeres?)
	zcat ${wd}/tumor_cov5.txt_f.gz | awk '{if ($8<200 && $9<200) {print $0}}' | gzip > ${wd}/tumor_cov5.txt_f_200std.gz 

	#split out the centromere clusters and other clusters, just for further analysis down the road
	zcat ${wd}/tumor_cov5.txt_f_200std.gz | awk '{if (NF>4) {print $0}}' | awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  | $b/bin/intersectBed -a - -b $cr | awk '{print $5"\t"$6"\t"$6+1"\t"$7"\t"$1"\t"$2"\t"$4"\t"$8"\t"$9"\t"$10}' | $b/bin/intersectBed -a - -b $cr  | awk '{print $5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$4"\t"$8"\t"$9"\t"$10}' | gzip > ${wd}/both_centro_tumor_cov5.txt_f.gz
	zcat ${wd}/tumor_cov5.txt_f_200std.gz | python $g/scripts/remove_lines.py ${wd}/both_centro_tumor_cov5.txt_f.gz |  gzip > ${wd}/not_both_centro_tumor_cov5.txt_f.gz
	#merge them back anyway
	zcat ${wd}/both_centro_tumor_cov5.txt_f.gz ${wd}/not_both_centro_tumor_cov5.txt_f.gz | gzip > ${wd}/all.txt_f.gz

	python $g/scripts/overlap.py 3000 ${wd}/normal_cov0.txt_f.gz ${wd}/all.txt_f.gz  |  awk '{if (NF>4) {print $0}}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}'  > $wd/nsubtract_centrosubtract_1000bp
	#python $g/scripts/overlap.py 3000 ${wd}/normal_cov0.txt_f.gz ${wd}/not_both_centro_tumor_cov5.txt_f.gz  |  awk '{if (NF>4) {print $0}}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}'  > $wd/nsubtract_centrosubtract_1000bp
	#python $g/scripts/overlap.py 50000 ${wd}/normal_cov0.txt_f.gz ${wd}/both_centro_tumor_cov5.txt_f.gz  |  awk '{if (NF>4) {print $0}}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}'  >> $wd/nsubtract_centrosubtract_1000bp

	#get the local bam files
	echo Getting local intervals...
	cat $wd/nsubtract_centrosubtract_1000bp | grep chr  | awk '{OFS="\t"; if ($3=="+") {print $1,$2,"-"} else {print $1,$2,"+"}; if ($6=="+") {print $4,$5,"+"} else {print $4,$5,"-"};}' | sort | uniq  | awk '{if ($3=="-") {print $1"\t"($2-1200)"\t"($2+200)} else {print $1"\t"($2-200)"\t"($2+1200)}}' | awk '{if ($2<=0) {$2=1} ; if ($3<=0) {$3=1}; print $1"\t"$2"\t"$3}' | $b/bin/sortBed -i - | $b/bin/mergeBed -i - > ${wd}/intervals

	function localBam {
		echo making local $ty bam
		local ty=$1
		$s view -H ${wd}/${ty}.bam > ${wd}/${ty}_sam_headers
		cat ${wd}/intervals | while read line; do
			echo $ty $line 1>&2
			$s view -H ${wd}/${ty}.bam | grep -i chr20 > /dev/null
			if [ $? -eq 0 ] ; then
				$s view ${wd}/${ty}.bam `echo $line | sed 's/\(chr\S*\)\s*\([0-9]*\)\s*\([0-9]*\)/\1:\2-\3/g'` 
			else
				$s view ${wd}/${ty}.bam `echo $line | sed 's/chr\(\S*\)\s*\([0-9]*\)\s*\([0-9]*\)/\1:\2-\3/g'` 
			fi
		done  | sort -S 5g | uniq | cat ${wd}/${ty}_sam_headers - | $s view -Sb -o ${wd}/${ty}_uniq.bam - 
	}


	for ty in $tys; do
		localBam $ty &
	done
	wait

	#refine the alignments and join them
	echo refining alignments...
	$s sort ${wd}/tumor_uniq.bam ${wd}/tumor_uniq_sorted
	tumor_mean=`cat ${wd}/tumor_mean_and_std.txt | awk '{print $1}'`
	tumor_std=`cat ${wd}/tumor_mean_and_std.txt | awk '{print $2}'`
	$s view ${wd}/tumor_uniq_sorted.bam | $g/clustering/cluster_em ${tumor_mean} ${tumor_std} ${wd}/nsubtract_centrosubtract_1000bp $ref > $wd/nsubtract_centrosubtract_1000bp_emd
	#cat $wd/nsubtract_centrosubtract_1000bp $wd/nsubtract_centrosubtract_1000bp_emd | sort | awk '{c=$1","$2","$3","$4","$5","$6; if (c==a) {s_new=$NF; if (s_new>s_old) {print c,s_new} else {print a,s_old}; a=""; c=""; } else {if (a!="") {print a,s_old}; a=c; s_old=$NF}}' | sed 's/,/\t/g' > $wd/nsubtract_centrosubtract_1000bp_joined

	cat $wd/nsubtract_centrosubtract_1000bp_emd | awk '{c=$1","$2","$3","$4","$5","$6; s=$NF; if (s>2) {print c,s}}' | sed 's/,/\t/g' > $wd/nsubtract_centrosubtract_1000bp_joined

	#now get the arc coverages
	echo getting arc coverages
	normal_mean=`cat ${wd}/normal_mean_and_std.txt | awk '{print $1}'`
	normal_std=`cat ${wd}/normal_mean_and_std.txt | awk '{print $2}'`
	cat $wd/nsubtract_centrosubtract_1000bp_joined | grep chr  | awk '{print $1"\t"$2; print $4"\t"$5; }' | sort | uniq > ${wd}/nsubtract_centrosubtract_1000bp_joined_intervals
	$s view $wd/normal_uniq.bam | $g/arc_coverage/arc_coverage ${normal_mean} ${normal_std} $wd/nsubtract_centrosubtract_1000bp_joined_intervals > $wd/nsubtract_centrosubtract_1000bp_joined_arc_coveage

	cat $wd/nsubtract_centrosubtract_1000bp_joined | sed 's/\(chr[^:]*\):\([0-9]*\)\([+-]\)\s\(chr[^:]*\):\([0-9]*\)\([+-]\)/\1\t\2\t\3\t\4\t\5\t\6/g' | awk '{OFS="\t"; type=0; if ($3=="+") {if ($6=="+") {type=0} else {type=2} } else { if ($6=="+") {type=3} else {type=1} }; print $1,$2,$5,type,$7,0,0,0.0,0,$4,"EDGE" }' > $wd/nsubtract_centrosubtract_1000bp_joined_links
else
	echo Skipping refinement of bps
fi


if [ ! -e ${wd}/hmm_out -o $wd/nsubtract_centrosubtract_1000bp_joined_links -nt ${wd}/hmm_out ] ; then
	echo running hmm
	$g/hmm_p $wd/nsubtract_centrosubtract_1000bp_joined_links ${wd}/tumor_cov ${wd}/normal_cov  > ${wd}/hmm_out
else
	echo skipping hmm run
fi

echo running walker
offsets="0 10 20 30"
multipliers="1 10 20"
flows="0 1 128" # 2 4 8 16 32 64 128"
for o in $offsets; do
for m in $multipliers; do
for flow in $flows; do
	d=$wd/flows/f${flow}/o${o}_m${m}/
	mkdir -p $d
	if [ 1 -lt 3 -o ! -e $d/walker_out -o ${wd}/hmm_out -nt $d/walker_out ]; then
		echo Forcing walker run
		rm $d/*
		pushd $d
		$g/walker/walker_pp $wd/nsubtract_centrosubtract_1000bp_joined $wd/nsubtract_centrosubtract_1000bp_joined_arc_coveage ${wd}/hmm_out $flow N $o $m > walker_out
		zcat problem_file | $c | gzip > solved.gz
		pypy $g/walker/flow_to_graph_v2.py problem_file.gz solved.gz `echo $m 2 | awk '{print $1/$2}'` $o 10000 | sed 's/,//g' > graph
		popd
	fi
done
done
done
