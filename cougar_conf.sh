g=$COUGARD
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/ # i dont think this is actually used anymore.... 
c=$g/dependencies/cs2/cs2.exe
cg=/usr/bin/cgquery
gt=/usr/bin/gtdownload
key=/dupa-filer/misko/tcga/cghub.key #required if you are downloading TCGA data
p=$g/dependencies/picard-tools-1.56/ #this should be in this repo.. 
gurobi=/dupa-filer/misko/gurobi/gurobi550/linux64/bin/gurobi_cl


#need this for gurobi
#export LD_LIBRARY_PATH=/home/buske/arch/sge6/lib:/home/buske/arch/sge6/lib::/dupa-filer/misko/gurobi/gurobi550/linux64/lib/:/dupa-filer/misko/gurobi/gurobi550/linux64/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/dupa-filer/misko/gurobi/gurobi550/linux64/lib/



if [ ! -d ${g} ] ; then
	echo "Failed to find GouGaR path installed... ${g} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi
if [ ! -f ${gurobi} ] ; then
	echo "Failed to find GUROBI installed... ${gurobi} is not valid"
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

