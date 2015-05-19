#!/usr/bin/python


import sys
import math
from copy import deepcopy

import gzip




node_to_pos={}
genomic_edges={}
somatic_edges={}



def to_pos(s):
	s=s.split(':')
	return (to_chr(s[0]),int(s[1]))

def to_chr(s):
	s=s.lower()
	if s[:3]=='chr':
		s=s[3:]
		if s=='x':
			return 23
		if s=='y':
			return 24
		if s=='m':
			return 25
	return int(s)

if len(sys.argv)!=4:
	print "%s problem_file contig_file contig" % sys.argv[0]
	sys.exit(1)


def read_problem_file(filename):
	h=gzip.open(filename,'r')
	for line in h:
		line=line.split()
		if line[1]=='NODE':
			#node line
			#c NODE 24 chr1:202709642
			n=int(line[2])
			node_to_pos[n]=to_pos(line[3])
		elif line[1]=='Somatic':
			#Somatic line
			#c Somatic       chr1:202869943  chr1:204056467  2   cap normal tumor
			f=to_pos(line[2])
			t=to_pos(line[3])
			if f>t:
				print >> sys.stderr, "ERROR" 
				sys.exit(1)
			typ=int(line[4])
			cost=int(line[5])
			cap=int(line[6])
			normal=int(line[7])
			tumor=int(line[8])
			if (f,t) not in somatic_edges:
				somatic_edges[(f,t)]=typ
			if (f,t) in genomic_edges:
				print >> sys.stderr, "FAILED PRECONDITION",line
				#sys.exit(1)
		elif line[1]=='Genomic':
			#genomid line
			#c Genomic       chr2:8917257    chr2:8921728 cap length normal tumor
			f=to_pos(line[2])
			t=to_pos(line[3])
			if f>t:
				print >> sys.stderr, "ERROR" 
				sys.exit(1)
			cost=int(line[4])
			cap=int(line[5])
			length=int(line[6])
			normal=int(line[7])
			tumor=int(line[8])
			if (f,t) not in genomic_edges:
				genomic_edges[(f,t)]=length
			if (f,t) in somatic_edges:
				print >> sys.stderr, "FAILED PRECONDITION"
				#sys.exit(1)
	h.close()
	

problem_filename=sys.argv[1]
simple_contigs_filename=sys.argv[2]
ci=int(sys.argv[3])

read_problem_file(problem_filename)


simple_contigs=[]
simple_contigs_f=open(simple_contigs_filename)
for line in simple_contigs_f:
        simple_contigs.append(eval(line.strip())[0])

c=simple_contigs[ci]
walk=[]
for x in range(len(c)):
	fn=c[x]
	tn=c[(x+1)%len(c)]
	fnp=node_to_pos[fn]
	tnp=node_to_pos[tn]
	if (fnp,tnp) not in genomic_edges and (tnp,fnp) not in genomic_edges:
		continue
	if len(walk)==0 or walk[-1][1]!=fnp:
		walk.append([fnp,tnp])
	else:
		walk[-1][1]=tnp

for (fn,tn) in walk:
	if (fnp,tnp) not in genomic_edges and (tnp,fnp) not in genomic_edges:
		#somatic edge
		print fn,tn,"*"
	else:
		print fn,tn	

