#60     /dupa-filer/misko/tcga/test/TCGA-ER-A19T/nsubtract_centrosubtract_1000bp_links  /dupa-filer/misko/tcga/test/TCGA-ER-A19T/tumor_cov      /dupa-filer/misko/tcga/test/TCGA-ER-A19T/normal_cov
#2       chr1:1  chr1:55157      55156   1588    857
#1       chr1:55157      chr1:98321      43164   1788    1524
#2       chr1:98321      chr1:246997     148676  3815    1921
#3       chr1:246997     chr1:539553     292556  3611    1305


import sys


l={}
c={}

for line in sys.stdin:
	if line[0]=='#':
		continue
	line=line.split()
	chr=line[1].split(':')[0]
	if chr not in l:
		l[chr]=0
	if chr not in c:
		c[chr]={}
	cn=int(line[0])
	if cn not in c[chr]:
		c[chr][cn]=0
	l[chr]+=int(line[3])
	c[chr][cn]+=int(line[3])


chrs=l.keys()
chrs.sort()


states=0
for chr in chrs:
	cns=c[chr].keys()
	cns.sort()
	xs=[]
	for cn in cns:
		x=float(c[chr][cn])/l[chr]
		if x>0.15:
			xs.append((x,str(cn)+":%0.2f" % x))
	xs.sort()
	#print "I"*len(xs)
	print chr,",".join(map(lambda x : x[1], xs))
	states+=len(xs)

print sys.argv[1],states
