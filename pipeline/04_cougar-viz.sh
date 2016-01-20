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

f=$1

if [ ! -d $f ] ; then
	echo "Folder does not exist! $f "
	exit
fi


if [ ! -f $f/ref ]; then
	echo "Cannot find file $f/ref , what is this hg18 or hg19?"
	exit
fi

ref=`cat $f/ref | head -n 1`

find $f | grep Qg_ | while read line; do
	if [ "${line##*.}" != "svg" ] ; then 
	o=$line".svg"
	if [ ! -f $o -o $line -nt $o ]; then 
#[misko@supa20 CouGaR]$ ../racket-6.2/bin//gracket cougarviz/cougar-viz.rkt 
#cougar-viz.rkt gene_annotations_file genomic-somatic_regions_file svg_out_filename
#[misko@supa20 CouGaR]$ ls cougarviz/gene_annotations/
#genes_wstrands_formatted_hg18.txt  genes_wstrands_formatted_hg19.txt
		$gracket $cougarviz/cougar-viz.rkt $cougarviz/gene_annotations/genes_wstrands_formatted_${ref}.txt $line $o
	fi
	fi
done



