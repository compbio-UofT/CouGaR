#/usr/bin/python

import sys
import gzip

#need to read in problem file for lower bounds on flow
#need to read in solution to simplify


edges_low={}

if len(sys.argv)!=3:
	print "%s problem_file.gz solved.gz" % (sys.argv[0])
	sys.exit(1)

h=gzip.open(sys.argv[1])
for line in h:
	if line[0]=='a':
		f,t,low,cap,score=map(int,line.split()[1:])
		if low!=0:
			if (f,t) not in edges_low:
				edges_low[(f,t)]=0
			edges_low[(f,t)]+=low


nodes=set()
edges={}

o=[]
o.append("c graph is simplified...")

h=gzip.open(sys.argv[2])
for line in h:
	if line[0]=='f':
		# a flow
		f,t,fl=map(int,line.split()[1:])
		if fl>0:
			#(f,t,low,cap,score)
			nodes.add(f)
			nodes.add(t)
			if (f,t) not in edges:
				edges[(f,t)]=0
			edges[(f,t)]+=fl


nnodes=0
nedges=0

for node in nodes:
	o.append("n\t%d\t0" % node)
	nnodes+=1	

for (f,t) in edges:
	low=0
	if (f,t) in edges_low:
		low=edges_low[(f,t)]
	if edges[(f,t)]>0:
		nedges+=1
		o.append("a\t%d\t%d\t%d\t%d\t%d" % (f,t,low,edges[(f,t)],-1))

print "p\tmin\t%d\t%d" % (max(nodes),nedges)
print "\n".join(o)

