#Chromosome      Start   End     Probe_Number    Segment_Mean
#1       554268  639581  3       -0.561
#1       736483  9985450 1500    0.5427
#1       9987720 10304808        95      0.2286
#1       10306679        10941525        115     0.5325
#1       10946936        15407023        651     -0.3389
#1       15409548        25457903        1816    -0.5529
#1       25482036        25529549        4       -1.6003
#1       25538212        27483586        409     -0.5603

#Chromosome	Start	End	Tumor_Count	Normal_Count	Segment_Mean
#chr1	10208	80595984	8866353	6918315	0.0244942024795493
#chr1	80595985	80670588	7362	7263	-0.391768643137613
#chr1	80670589	249240606	16753943	13041317	-0.00408394286779093
#chr2	10001	78614929	8674657	7462512	-0.140789416365168
#chr2	78614930	136507121	5575636	4348300	0.00181119735833432
#chr2	136507122	136670882	23904	15143	0.31581618660609

import sys

if len(sys.argv)!=5:
	print sys.argv[0],"fname length parts spike"
	sys.exit(2)

fn=sys.argv[1]
l=int(sys.argv[2])
p=int(sys.argv[3])
s=float(sys.argv[4])


segp=0

h=open(sys.argv[1])
for line in h:
	if line[0]=='C' or line[0]=='S':
		continue
	line=line.strip().split()
	ln=len(line)
	segl=int(line[2])-int(line[1])
	segs=float(line[ln-1])
	if segl>l and segs>s:
		segp+=1
h.close()

if segp<p:
	print "F",sys.argv[1]
	sys.exit(1)
else:
	print "T",sys.argv[1]
	sys.exit(0)
