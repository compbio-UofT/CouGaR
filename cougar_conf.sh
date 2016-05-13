#!/bin/bash

if [ -z "${g}" ] ; then
	echo "Failed to find GouGaR path installed... ${g} is not valid"
	echo "Trying current folder " `pwd`
        COUGARD=`pwd`/
fi

g=$COUGARD

#USER NEEDS TO INSTALL THIS AND FILL IN THIS PATH
#gurobi=/dupa-filer/misko/gurobi/gurobi650/linux64/bin/gurobi_cl
gurobi_install_path=/home/misko/gurobi651/
gurobi=${gurobi_install_path}/linux64/bin/gurobi_cl
j=/usr/bin/java
#need this for gurobi
#export LD_LIBRARY_PATH=/home/buske/arch/sge6/lib:/home/buske/arch/sge6/lib::/dupa-filer/misko/gurobi/gurobi550/linux64/lib/:/dupa-filer/misko/gurobi/gurobi550/linux64/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${gurobi_install_path}/linux64/lib/

#needed if you are trying to run COUGARVIZ
gracket=`which gracket`


#USED ONLY IF YOU ARE TRYING TO DOWNLOAD TCGA CONTENT WITH YOUR OWN KEYS
#cg=/usr/bin/cgquery
#gt=/usr/bin/gtdownload
#key=/dupa-filer/misko/tcga/cghub.key #required if you are downloading TCGA data

s=`which samtools`
#b=/filer/misko/bedtools-2.17.0/ # i dont think this is actually used anymore.... 
c=$g/dependencies/cs2/cs2.exe
p=$g/dependencies/picard-tools-1.56/ #this should be in this repo.. 

#programs that can be autosetup if they are missing
hgref=$g/refs
cougarviz=$g/cougarviz



if [ ! -f "${g}/cougar_conf.sh" ] ; then
	echo "Failed to find GouGaR path installed... ${g} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
fi


if [ ! -d "${cougarviz}" ] ; then
	echo "Failed to find GouGaR-viz path installed... ${cougarviz} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
        echo "downloading cougar-viz in 5 seconds..."
        sleep 5
	pushd $g
	git clone https://github.com/compbio-UofT/CouGaR-viz.git cougarviz
        popd
	if [ ! -d "${cougarviz}" ] ; then
		echo "Failed to git clone cougarviz repo.. mising git ? or missing internet?"
		exit
	fi
fi
if [ -d "${cougarviz}" ] ; then
	echo "Found cougar-viz"
fi


## FIND GUROBI
if [ ! -f ${gurobi} ] ; then
	echo "Failed to find GUROBI installed... ${gurobi} is not valid"
	echo "Please install Gurobi and fill in the path in this file...."
	echo "http://www.gurobi.com/downloads/download-center"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
else 
	echo "Found gurobi"
fi

#test if gurobi works
${gurobi} ResultFile=${g}/cplex/test.sol ${g}/cplex/test.lp  > ${g}/gurobi.log
diff ${g}/cplex/test.sol ${g}/cplex/test.sol
if ! cmp -s "${g}/cplex/test.sol" "${g}/cplex/test.sol_correct"
then
	echo "Gurboi does not seem to be working.. please check ${g}/gurobi.log for more details"
	echo "Maybe the library path has not been set for gurobi? please set this in cougar_conf.sh"
	#exit
else 
	echo "Passed gurobi test"
fi

#see if racket is installed
if [ ! -f ${gracket} ] ; then
	echo "Failed to find RACKET installed... ${gracket} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
else 
	echo "Found gracket"
fi


#see if samtools is installed, if not try to compile it
if [ -z "${s}" -o ! -f "${s}" ] ; then
	if [ ! -f "${g}/samtools-1.3.1/samtools" ] ; then
		echo "Failed to find SAMTOOLS installed... ${s} is not valid"
		echo "Please update $COUGARD/cougar_conf.sh"
		echo "trying to compile samtools in 5 seconds..."
		pushd $g
        	wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
        	bunzip2 samtools-1.3.1.tar.bz2
        	tar -xvf samtools-1.3.1.tar
        	pushd samtools-1.3.1
		make
		popd
	fi
	s=${g}/samtools-1.3.1/samtools
	if [ ! -f "${s}" ] ; then
		echo "seriously failed to find samtools or compile it..."
		exit
	fi	
fi

if [  -f "${s}" ] ; then
	echo "Found samtools"
fi


#check if java is installed
if [ ! -f ${java} ] ; then
	echo "Failed to find SAMTOOLS installed... ${java} is not valid"
	echo "Please update $COUGARD/cougar_conf.sh"
	exit
else
	echo "Found java"
fi

#check for CS2 installation
if [ ! -f ${c} ] ; then
	pushd $COUGARD/dependencies/cs2 
	make
	popd
	if [ ! -f ${c} ] ; then
		echo "Failed to find CS2 installed... ${c} is not valid"
		echo "please check the folder $COUGARD/dependencies/cs2 and run make"
		echo "Please update $COUGARD/cougar_conf.sh"
		exit
	fi
fi
if [ -f ${c} ] ; then
	echo "Found CS2"
fi



hgref=$g/refs


md5sum --check $g/refs/checksums.txt
if [ $? -ne 0 ]; then
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
md5sum --check $g/refs/checksums.txt
if [ $? -ne 0 ]; then
	echo "Failed to download or find reference fasta files..."
	exit
else
	echo "Found HG18 and HG19 reference files!"
fi
