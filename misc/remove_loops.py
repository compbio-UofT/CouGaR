#!/usr/bin/python

import os
import sys
from copy import deepcopy
import random
#from multiprocessing import Pool


from subprocess import Popen, PIPE, STDOUT


import gzip
import time

params={'num_nodes':0,'k':0,'cost':0,'mins':0}
edges={}

full_params={'num_nodes':0,'k':0,'cost':0,'mins':0}
full_edges={}

pid=os.getpid()

#provide extra edge to location from node 3
def print_graph(flow_restriction={},start_node=None,extra=-1):
	sx=0
	arc_lines=[]
	for f in edges:
		for t in edges[f]:
			score,low,cap=edges[f][t]
			used=0
			if f in flow_restriction and t in flow_restriction[f]:
				used=flow_restriction[f][t]
			if used!=cap:
				low=max(0,low-used)
				arc_lines.append( "a\t%d\t%d\t%d\t%d\t%d" % (f,t,low,cap-used,score))
				sx+=-1*(cap-used)
	if extra>0:
		arc_lines.append( "a\t%d\t%d\t%d\t%d\t%d" % (start_node,extra,1,1,0))
	arc_lines.append("a\t%d\t%d\t%d\t%d\t%d" % (1,1,0,0,0))

	lines=[]
	lines.append("p\tmin\t%d\t%d" % (params['num_nodes'],len(arc_lines)))
	for x in range(params['num_nodes']):
		lines.append("n\t%d\t0" % (x+1))
	for arc_line in arc_lines:
		lines.append(arc_line)
	return sx,"\n".join(lines)

def get_flow(s):
	fld={}
	cost=1337
	#stream = os.popen("cat " + problem_filename + " | /data/misko/2013.04.12/cs2-4.3/cs2.exe")
	
	p = Popen("/filer/misko/cs2.exe",stdin=PIPE,stdout=PIPE)
	stream=p.communicate(input=s)[0].splitlines()
	for line in stream:
		if line.find('cost')>0:
			cost = int(line.split()[4])
		if line[0]=='f':
			line=line.split()
			f=int(line[1])
			t=int(line[2])
			fl=int(line[3])
			if fl>0:
				if f not in fld:
					fld[f]={}
				if t not in fld[f]:
					fld[f][t]=0
				fld[f][t]+=fl
	#print >> sys.stderr, " SOLVED FLOW ", cost
	return fld,cost

def read_graph(x,l_edges,l_params,full=False):
	#expect graph
	walks=0
	for line in x:
		if line[0]=='c':
			#comment, essentially skip
			pass
		elif line[0]=='n':
			#handle new node
			pass
		elif line[0]=='a':
			#handle a arc
			line=line.split()
			f=int(line[1])
			t=int(line[2])
			low=int(line[3])
			walks=max(walks,low)
			cap=int(line[4])
			score=int(line[5])
			if f not in l_edges:
				l_edges[f]={}
			if t not in l_edges[f]:
				l_edges[f][t]=[]
			elif not full:
				print "already read this edge"
				sys.exit(1)
			if full:
				l_edges[f][t].append((score,low,cap))
			else:
				l_edges[f][t]=(score,low,cap)
		elif line[0]=='p':
			#handle program parameters
			line=line.split()
			l_params['num_nodes']=int(line[2])
			l_params['num_arcs']=int(line[3])
	l_params['walks']=walks
	return l_edges,l_params




if len(sys.argv)!=3:
	print "%s full_pb.gz simple_pb.gz" % sys.argv[0]
	sys.exit(1)

fpb_fname=sys.argv[1]
spb_fname=sys.argv[2]

print "Reading in simple graph...",
read_graph(gzip.open(spb_fname),edges,params)
print params
print "Reading in full graph...",
read_graph(gzip.open(fpb_fname),full_edges,full_params,True)
print full_params


#reduce the scores that we dont need
keep={}
for f in full_edges:
	keep[f]={}
	for t in full_edges[f]:
		full_edges[f][t].sort()
		keep[f][t]=0
		if f in edges and t in edges[f]:
			keep[f][t]=edges[f][t][2]

for f in keep:
	for t in keep[f]:
		full_edges[f][t]=full_edges[f][t][:keep[f][t]]
		
				

def cloop(l):
	return tuple(l[l.index(min(l)):]+l[:l.index(min(l))])

def dfs(flr,found):
	nodes=set()
	for f in edges:
		nodes.add(f)
		for t in edges[f]:
			nodes.add(t)
	rnodes=list(nodes)
	random.shuffle(rnodes)
	pi={}
	seen=set()
	for v in rnodes:
		q=[([],v)]
		while len(q)>0:
			parents,node=q.pop()
			try:
				for t in edges[node]:
					used=0
					if node in flr and t in flr[node]:
						used=flr[node][t]
					if edges[node][t][2]>used:
						if t in parents:
							#cycle?
							p=cloop(parents[parents.index(t):]+[node])
							if p not in found:
								return p
							else:
								#ignore edge
								pass
						elif t not in seen:
							seen.add(t)
							q.append((parents+[node],t))
			except:
				print parents,node
				sys.exit(1)
	return None



found={}
flr={}
for f in edges:
	flr[f]={}
	for t in edges[f]:
		flr[f][t]=0
k=dfs(flr,found)
while k!=None:
	#try to check for sub loop
	for x in k:
		if k.count(x)>1:
			print "LOOPY"
			sys.exit(1)
	#try to make new edge list
	tk=tuple(k)
	if tk not in found:
		found[tk]=0
	else:
		print "EERROROROROROOROROROROR"
		sys.exit(1)
	cap_ok=True
	m=0
	while cap_ok:
		flr_p=deepcopy(flr)
		for x in range(len(k)):
			f=k[x]
			t=k[(x+1)%len(k)]
			#make sure that we can subtract this based on cap
			used=flr_p[f][t]
			if edges[f][t][2]>used:
				flr_p[f][t]+=1
			else:
				cap_ok=False
				break
		if cap_ok:
			#lets try to flow solve
			#print "trying to solve"
			max_cost,g=print_graph(flr_p)
			local_fl,flow_cost=get_flow(g)
			if flow_cost==max_cost:
				m+=1
				#print "upping cycle count",m
				flr.update(flr_p)
			else:
				#print "cannot up cycles any more",m
				cap_ok=False
		else:
			#print "no good cap"
			max_cost,g=print_graph(flr)
			local_fl,flow_cost=get_flow(g)
			#print max_cost,flow_cost
		found[tk]=m
	k=dfs(flr,found)

def score(k,m):
	score=0
	for x in range(len(k)):
		f=k[x]
		t=k[(x+1)%len(k)]
		score+=sum(map( lambda z : z[0], full_edges[f][t][-m:]))
 	return score	

outs=[]
for k in found:
	s=score(k,found[k])
	x=""
	if 1 in k or 0 in k:
		x="*"
	outs.append((s,len(k),found[k],x))
outs.sort()
for out in outs:
	print "\t".join(map(str,out))
max_cost,g=print_graph(flr)
local_fl,flow_cost=get_flow(g)
print max_cost, flow_cost, local_fl
