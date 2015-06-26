g=/dupa-filer/misko/cougar/CouGaR/
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/
c=/data/misko/2013.04.12/cs2-4.3/cs2.exe
cg=/usr/bin/cgquery
gt=/usr/bin/gtdownload
key=/dupa-filer/misko/tcga/cghub.key
p=/filer/misko/picard/picard-tools-1.56/
gurobi=/dupa-filer/misko/gurobi/gurobi550/linux64/bin/gurobi_cl

#need this for gurobi
#export LD_LIBRARY_PATH=/home/buske/arch/sge6/lib:/home/buske/arch/sge6/lib::/dupa-filer/misko/gurobi/gurobi550/linux64/lib/:/dupa-filer/misko/gurobi/gurobi550/linux64/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/dupa-filer/misko/gurobi/gurobi550/linux64/lib/
