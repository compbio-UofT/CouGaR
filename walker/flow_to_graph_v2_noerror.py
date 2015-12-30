#!/usr/bin/python


import sys
import math
from copy import deepcopy

import gzip




#WARNING MAKING ASSUMPTIONS ABOUT EDGES.......
#BAD ASSUMPTIONS


keep_empty=False # set this to true to keep edges with no flow!

node_to_pos={}
genomic_edges={}
somatic_edges={}

expected={}

break_points=set()

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
			#print >> sys.stderr, "chrM is enabled!"
			return 25
	return int(s)

if len(sys.argv)!=6:
	print "%s problem_file flow scale offset smooth" % sys.argv[0]
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
			if (f,t) not in expected:
				expected[(f,t)]=(normal,tumor)
			break_points.add(f)
			break_points.add(t)
			if (f,t) in genomic_edges:
				print >> sys.stderr, "FAILED PRECONDITION, TODO fix this",line
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
			if (f,t) not in expected:
				expected[(f,t)]=(normal,tumor)
			if (f,t) not in genomic_edges:
				genomic_edges[(f,t)]=length
			if (f,t) in somatic_edges:
				print >> sys.stderr, "FAILED PRECONDITION, TODO fix this"
				#sys.exit(1)
	h.close()
	
def annotate_flow(p):	
	used={}
	for fn,tn,l in p:
		if fn>2 and tn>2:
			f=node_to_pos[fn]
			t=node_to_pos[tn]
			# not the source and sink chromosomes
			if f[0]>0 and t[0]>0:
				#get the canonical
				c=(f,t)
				if f>t:
					c=(t,f)
				if c[0]>c[1]:
					print >> sys.stderr, "FAIL"
					sys.exit(1)
				if c not in used:
					used[c]={'genomic':0,'somatic':0}
				#find out if we are using somatic or genomic
				#if l>1:
				#	print fn,tn,l,f,t,c
				for x in range(l):
					if c in genomic_edges:
						used[c]['genomic']+=1
					else:
						used[c]['somatic']+=1
	return used



def read_flow_file(filename):
	#f       1       3         45
	#f       3      65          0
	#f       3     117          0
	h=gzip.open(filename,'r')
	flows=[]
	for line in h:
		line=line.split()
		if line[0]=='f':
			#its a flow line, lets assign flow
			f=int(line[1])
			t=int(line[2])
			fl=int(line[3])
			flows.append((f,t,fl))					
	h.close()
	print "; " + filename 
	used=annotate_flow(flows)
	
	used_breakpoints=set()

	l_somatic_edges=[]
	for f,t in somatic_edges:
		u=0
		if (f,t) in used:
			u=used[(f,t)]['somatic']*scale
			if u>0:
				u+=offset
			if abs(u-int(u))<0.00001:
				u=int(u)
			if used[(f,t)]['somatic']>0 and u<=0:
				print >> sys.stderr, "ERROR 3235",used[(f,t)]['somatic']
				sys.exit(1)
			if keep_empty or used[(f,t)]['somatic']>0.01:
				l_somatic_edges.append((f,t,somatic_edges[(f,t)],u,expected[(f,t)][0],expected[(f,t)][1]))
				used_breakpoints.add(f)
				used_breakpoints.add(t)
		#if (t[0]=='chr21' or t[0]=='21' or t[0]==21) and 14366000<t[1]:
		#	print >> sys.stderr, (f,t), u, used[(f,t)]['somatic'],'SOM'

	l_genomic_edges=[]
	genomic_edges_keys=genomic_edges.keys()
	genomic_edges_keys.sort()
	for f,t in genomic_edges_keys:
		u=0
		if (f,t) in used:
			u=used[(f,t)]['genomic']*scale
			if u>0:
				u+=offset
			if abs(u-int(u))<0.00001:
				u=int(u)
		if f not in used_breakpoints and len(l_genomic_edges)>0 and l_genomic_edges[-1][1]==f and l_genomic_edges[-1][3]==u:
			#add it to the previous
			prev=l_genomic_edges[-1]
			l_genomic_edges[-1]=(prev[0],t,0,u,prev[4]+expected[(f,t)][0],prev[5]+expected[(f,t)][1])
		else:
			if len(l_genomic_edges)>0 and l_genomic_edges[-1][0][0]==f[0] and 0<f[1]-l_genomic_edges[-1][1][1]<smooth:
				l_genomic_edges.append((l_genomic_edges[-1][1],f,0,0,-1,-1))
			l_genomic_edges.append((f,t,0,u,expected[(f,t)][0],expected[(f,t)][1]))
		#if (f[0]=='chr21' or f[0]=='21' or f[0]==21) and 14366000<f[1]:
		#	print >> sys.stderr, (f,t), u, used[(f,t)]['somatic']


	sz=len(l_genomic_edges)+1
	while sz>len(l_genomic_edges):
		print >> sys.stderr, "ITERATION"
		sz=len(l_genomic_edges)
		if not keep_empty and len(l_genomic_edges)>1:
			to_remove=[]
			#check front
			(f,t,y,u,no,tu)=l_genomic_edges[0]
			(nf,nt,ny,nu,pno,ptu)=l_genomic_edges[1]
			if u==0 and f not in used_breakpoints and t not in used_breakpoints and t!=nf:
				to_remove.append((f,t,y,u,no,tu))
			#check back
			(f,t,y,u,no,tu)=l_genomic_edges[-1]
			(pf,pt,py,pu,pno,ptu)=l_genomic_edges[1]
			if u==0 and f not in used_breakpoints and t not in used_breakpoints and f!=pt:
				to_remove.append((f,t,y,u,no,tu))
			#check middle
			for x in range(1,len(l_genomic_edges)-1):
				(pf,pt,py,pu,pno,ptu)=l_genomic_edges[x-1]
				(f,t,y,u,no,tu)=l_genomic_edges[x]
				(nf,nt,ny,nu,nno,ntu)=l_genomic_edges[x+1]
				#if (f[0]=='chr10' or f[0]=='10' or f[0]==10) and 42529526==f[1]:
				#	print >> sys.stderr, "MATCH", l_genomic_edges[x-1],l_genomic_edges[x],l_genomic_edges[x+1]
				if u==0 and f not in used_breakpoints and t not in used_breakpoints and (f!=pt or t!=nf):
					to_remove.append((f,t,y,u,no,tu))
				#elif u==0 and nf!=t and t not in used_breakpoints:
				#	to_remove.append((f,t,y,u,e))
				#elif u==0 and f!=pt and f not in used_breakpoints:
				#	to_remove.append((f,t,y,u,e))
			for x in to_remove:
				l_genomic_edges.remove(x)
			
	print str(l_genomic_edges).replace(',','')
	print str(l_somatic_edges).replace(',','')

problem_filename=sys.argv[1]
#paths_filename=sys.argv[2]
flow_filename=sys.argv[2]
scale=float(sys.argv[3])
offset=float(sys.argv[4])
smooth=int(sys.argv[5])


read_problem_file(problem_filename)
#read_paths_file(paths_filename)
read_flow_file(flow_filename)

