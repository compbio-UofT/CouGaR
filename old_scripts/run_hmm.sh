#!/bin/bash

g=/filer/misko/mini_chr/git/minichr/

b=$1
d=`dirname ${b}`
f=`echo ${1} | sed "s/.*\/\([^/]*\)/\1/g" | sed 's/.bam//g'`


$g/hmm ${d}/${f}_clusters/*q15*cov3*links ${d}/tumor_cov ${d}/normal_cov > ${1}_hmm_out 
