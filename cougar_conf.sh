g=$COUGARD
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/ # i dont think this is actually used anymore.... 
c=$g/dependencies/cs2/cs2.exe
cg=/usr/bin/cgquery
gt=/usr/bin/gtdownload
key=/dupa-filer/misko/tcga/cghub.key #required if you are downloading TCGA data
p=$g/dependencies/picard-tools-1.56/ #this should be in this repo.. 
gurobi=/dupa-filer/misko/gurobi/gurobi650/linux64/bin/gurobi_cl
gracket=/dupa-filer/misko/tcga/paper/gits/racket-6.2/bin/gracket
cougarviz=$g/cougarviz

hgref=$g/refs



#need this for gurobi
#export LD_LIBRARY_PATH=/home/buske/arch/sge6/lib:/home/buske/arch/sge6/lib::/dupa-filer/misko/gurobi/gurobi550/linux64/lib/:/dupa-filer/misko/gurobi/gurobi550/linux64/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/dupa-filer/misko/gurobi/gurobi650/linux64/lib/




if [ ! -d "${g}" ] ; then
	echo "Failed to find GouGaR path installed... ${g} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi
if [ ! -d "${cougarviz}" ] ; then
	echo "Failed to find GouGaR-viz path installed... ${cougarviz} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi
if [ ! -f ${gurobi} ] ; then
	echo "Failed to find GUROBI installed... ${gurobi} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi
#test if gurobi works
${gurobi} ResultFile=${g}/cplex/test.sol ${g}/cplex/test.lp  > ${g}/gurobi.log
diff ${g}/cplex/test.sol ${g}/cplex/test.sol
if ! cmp -s "${g}/cplex/test.sol" "${g}/cplex/test.sol_correct"
then
	echo "Gurboi does not seem to be working.. please check ${g}/gurobi.log for more details"
	echo "Maybe the library path has not been set for gurobi? please set this in cougar_conf.sh"
	exit
fi
if [ ! -f ${gracket} ] ; then
	echo "Failed to find RACKET installed... ${gracket} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi
if [ ! -f ${s} ] ; then
	echo "Failed to find SAMTOOLS installed... ${s} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi
if [ ! -f ${java} ] ; then
	echo "Failed to find SAMTOOLS installed... ${java} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi
if [ ! -f ${c} ] ; then
	echo "Failed to find CS2 installed... ${c} is not valid"
	echo "please check the folder $COUGARD/dependencies/cs2 and run make"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi



hgref=$g/refs

if [ ! -d $hgref ]; then
	mkdir $hgref
	pushd $hgref
	chrs="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M"
	hgs="18 19"
	for hg in $hgs ; do 
		for chr in $chrs ; do 
			wget -O - http://hgdownload.cse.ucsc.edu/goldenpath/hg${hg}/chromosomes/chr${chr}.fa.gz | gunzip 
		done | gzip > hg${hg}.fa.gz
	done
	popd
fi
