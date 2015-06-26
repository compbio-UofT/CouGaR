#!/bin/bash


if [ $# -ne 1 ]; then
	echo $0 downloadable.xml
	exit
fi 

f=$1

cat $1 | grep disea | head -n 1 | sed 's@.*<disease_abbr>\(.*\)</disease_abbr>.*@\1@g'
