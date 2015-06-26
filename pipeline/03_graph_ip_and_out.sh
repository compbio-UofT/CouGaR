#!/bin/bash

#check that COUGARD has been properly set
if [ -z "$COUGARD" ]; then echo "Please set COUGAR directory COUGARD=path"; exit; fi
if [ ! -d "$COUGARD" ] ; then echo "COUGARD enviornment variable is not set properly"; exit; fi

#load the configuation 
. $COUGARD/cougar_conf.sh

if [ $# -ne 2 ]; then 
	echo $0 folder tcga-id
	exit
fi
wd=$1
id=$2

function mwalker {
	i=$1
	sq=$2
	m=$3
	pushd $wd
	#sed out the chr25
	sed -i 's/chr25/chrM/g' hmm
	sed -i 's/chr25/chrM/g' nsub*

	#remake the mqs	 - just in case they werent right from before?
        #ref=`python $g/get_ref.py downloadable.xml | awk '{print $NF}' | head -n 1`
	ref=`cat ref`
        $g/mapability/bigWigAverageOverBed $g/mapability/$ref/* nsubtract_centrosubtract_1000bp.bed nsubtract_centrosubtract_1000bp.bed.out
        cat nsubtract_centrosubtract_1000bp.bed.out | awk '{print $1,$NF}' > nsubtract_centrosubtract_1000bp_mqs 

	$g/walker/walker_new $wd/nsubtract_centrosubtract_1000bp ${wd}/hmm N 0 ${m} 2 > walker_out_q${i}sq${sq}_m${m}
	mv problem_file.gz problem_file_q${i}sq${sq}_m${m}.gz
	zcat  problem_file_q${i}sq${sq}_m${m}.gz | /data/misko/2013.04.12/cs2-4.3/cs2.exe | gzip > solved_q${i}sq${sq}_m${m}.gz
	
	##########################
	## FLOW   PROBLEM   ######
	##########################
	#make the lp problem
	pypy $g/cplex/graph_c_to_lp.py problem_file_q${i}sq${sq}_m${m}.gz el_q${i}sq${sq}_m${m}.gz nsubtract_centrosubtract_1000bp_mqs ${i} ${sq} ${m} Qproblem_file_q${i}sq${sq}_m${m} > lp_prob_q${i}sq${sq}_m${m}.lp
	rm Qproblem_file_q${i}sq${sq}_m${m}.gz
	gzip Qproblem_file_q${i}sq${sq}_m${m}
	#solve the flow problem - using cs2.exe 
	zcat Qproblem_file_q${i}sq${sq}_m${m}.gz | $c > Qproblem_file_q${i}sq${sq}_m${m}.solved
	# make sure we remove if this has been run a second time
	if [ -e Qproblem_file_q${i}sq${sq}_m${m}.solved.gz ] ; then 
		rm  Qproblem_file_q${i}sq${sq}_m${m}.solved.gz
	fi
	gzip  Qproblem_file_q${i}sq${sq}_m${m}.solved
	#flow to graph formath - this just makes the solution human readable and plottable by cougar-viz
	pypy $g/walker/flow_to_graph_v2.py Qproblem_file_q${i}sq${sq}_m${m}.gz  Qproblem_file_q${i}sq${sq}_m${m}.solved.gz ${m} 0 1000 > Qg_lp_q${i}sq${sq}_m${m}
	#update graph with no flows - insert flows from the hmm solution to see how the hmm solution compares
	pypy $g/walker/insert_no_flows.py Qg_lp_q${i}sq${sq}_m${m} hmm > Qg_lp_q${i}sq${sq}_m${m}_whmm

	##########################
	## IP   PROBLEM   ########
	##########################
	#docompose the graph into contigs then, create and run the ip solver
	zcat Qproblem_file_q${i}sq${sq}_m${m}.solved.gz  > Qproblem_file_q${i}sq${sq}_m${m}.solved.tmp
	pypy $g/decompose/decompose_and_report.py Qproblem_file_q${i}sq${sq}_m${m}.solved.tmp decompositions_q${i}sq${sq}_m${m}.loops decompositions_q${i}sq${sq}_m${m}.loops.full
	pypy $g/cplex/graph_c_to_lp_wsimple.py problem_file_q${i}sq${sq}_m${m}.gz el_q${i}sq${sq}_m${m}.gz nsubtract_centrosubtract_1000bp_mqs ${i} 20 ${m} Qproblem_file_q${i}sq${sq}_m${m}.tmp decompositions_q${i}sq${sq}_m${m}.loops > ip_prob_q${i}sq${sq}_m${m}.lp

	#wait for IP solver - used when only 1 license available on another machine that montiros folder for '*.go' files
	#touch ip_prob_q${i}sq${sq}_m${m}.lp.go
	#while [ -e ip_prob_q${i}sq${sq}_m${m}.lp.go ] ; do
	#	sleep $[ ( $RANDOM % 25 )  + 1 ]
	#done

	#used to run Gurobi locally
	${gurobi} MIPGap=0 ResultFile=ip_prob_q${i}sq${sq}_m${m}.lp.sol ip_prob_q${i}sq${sq}_m${m}.lp
	cat ip_prob_q${i}sq${sq}_m${m}.lp.sol  | grep c | awk '{if ($2>0) {print $0}}'  | pypy $g/cplex/contigs_to_flow.py decompositions_q${i}sq${sq}_m${m}.loops | gzip > ip_prob_q${i}sq${sq}_m${m}.lp.flow.gz
	pypy $g/walker/flow_to_graph_v2_noerror.py Qproblem_file_q${i}sq${sq}_m${m}.gz ip_prob_q${i}sq${sq}_m${m}.lp.flow.gz ${m} 0 1000 > Qg_ip_q${i}sq${sq}_m${m}

	#generate a graph file that is human readable and usable by cougar-viz from the IP solution 
	cat ip_prob_q${i}sq${sq}_m${m}.lp.sol  | awk '{if ($2>0) {print $0}}' | grep "^c" | while read line; do 
		z=`echo $line | awk '{print $1}'`_m`echo $line | awk '{print $2}'`
		echo $line | pypy $g/cplex/contigs_to_flow.py decompositions_q${i}sq${sq}_m${m}.loops | gzip > q${i}sq${sq}_m${m}_${z}.gz
		pypy $g/walker/flow_to_graph_v2_noerror.py Qproblem_file_q${i}sq${sq}_m${m}.gz q${i}sq${sq}_m${m}_${z}.gz ${m} 0 1000 > g_q${i}sq${sq}_m${m}_${z}
	done

	rm Qproblem_file_q${i}sq${sq}_m${m}.solved.tmp

	##########################
	## STATS          ########
	##########################
	#ref=`python $g/get_ref.py  downloadable.xml  | awk '{print $NF}' | head -n 1`
	ref=`cat ref`
	#group=`$g/get_group.sh downloadable.xml`
	group=`cat subset`
	echo $id $group > group
	/usr/bin/pypy $g/walker/contigs_to_coords.py Qproblem_file_q${i}sq${sq}_m${m}.gz decompositions_q${i}sq${sq}_m${m}.loops.full ip_prob_q${i}sq${sq}_m${m}.lp.sol ${m} $g/genes/genes_wstrands_formatted_${ref}.txt $g/genes/oncos ${id} group > stats

	popd 
}

#run overlaps
#mwalker 4000 1000

#mwalker 1350 300 3 # contigQ somatic_base_Q multiplier 

#this was run with ip SQ = SQ / 2
#mwalker 1350 290 3 # contigQ somatic_base_Q multiplier 

#this was run with ip SQ = -20
#mwalker 1350 291 3 # contigQ somatic_base_Q multiplier 

#this was run with ip SQ = 20
mwalker 1350 292 3 # contigQ somatic_base_Q multiplier 
