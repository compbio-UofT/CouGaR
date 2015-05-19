#/usr/bin/python


import sys


if len(sys.argv)!=3:
	print "%s f1 f2 - result is f1-f2" % sys.argv[0]
	sys.exit(1)


f1fn=sys.argv[1]
f2fn=sys.argv[2]


def read_flow(fn):
	d={}
	f=open(fn,'r')
	for line in f:
		if line[0]=='f':
			line=line.split()
			s=int(line[1])
			t=int(line[2])
			l=int(line[3])
			if not s in d:
				d[s]={}
			if t not in d[s]:
				d[s][t]=0
			d[s][t]+=l		
	f.close()
	return d


f1=read_flow(f1fn)
f2=read_flow(f2fn)


for s in f2:
	if not s in f1:
		print >> sys.stderr, "missing ", s
	else:
		for t in f2[s]:
			if t not in f1[s]:
				print >> sys.stderr, "missing ", s,t
			else:
				f1[s][t]-=f2[s][t]
er=False
for s in f1:
	for t in f1[s]:
		if f1[s][t]<0 and not er:
			er=True
			print >> sys.stderr, "NEGATIVE VALUE!"
		f1[s][t]=max(0,f1[s][t])
		if f1[s][t]>=0:
			print "f\t%d\t%d\t%d" % (s,t,f1[s][t])
