#!/usr/bin/python

import os
import sys
from copy import deepcopy
import random
#from multiprocessing import Pool
from subprocess import Popen, PIPE, STDOUT

import time

start_node=0
params={'num_nodes':0,'k':0,'cost':0,'mins':0}
edges={}

pid=os.getpid()

#provide extra edge to location from node 3
def print_graph(flow_restriction={},extra=-1):
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

def read_graph(x):
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
			if f not in edges:
				edges[f]={}
			if t not in edges[f]:
				edges[f][t]=[]
			else:
				"already read this edge"
				sys.exit(1)
			edges[f][t]=(score,low,cap)
		elif line[0]=='p':
			#handle program parameters
			line=line.split()
			params['num_nodes']=int(line[2])
			params['num_arcs']=int(line[3])
	return walks

def get_node(cc, flow):
	nodes=rk[cc].keys()
	random.shuffle(nodes)
	for node in nodes:
		for x in edges[node]:
			score,low,cap=edges[node][x]
			used=flow[node][x]
			if used<cap:
				return node
	return None	

def get_candidates(last_node,flow):
	r=[]
	#want to return edges leaving the last node that have positive flow
	for node in edges[last_node]:
		score,low,cap=edges[last_node][node]
		used=0
		if last_node in flow and node in flow[last_node]:
			used=flow[last_node][node]
		if used<cap:
			r.append(node)
	random.shuffle(r)
	return r



max_so_far=[-100000000000]
go_back=0

def search(t):
	l,fl,path,to_use=t
	if to_use<3000 and to_use>-3000:
		print t
	#fname='/dev/shm/out_t'+str(os.getpid())
	fname='/dev/shm/out_%d' % pid
	#ouc=print_graph(fname+"Y",fl,extra=path[-1])
	#ofl,ocost=get_flow(fname+"Y")
	
	results=[]
	#keep going until we dont hit the sink!
	last_node = path[-1]
	#if (last_node==5 or last_node==6) and params['walks']==(path.count(5)+path.count(6)):
	#local_fl,cost=get_flow(fname)
	if last_node==start_node and to_use==0:
		unused_cost,g=print_graph(fl,extra=last_node)
		print >> sys.stderr, "DONE!",params['walks'],params['mins'],params['cost'],to_use,unused_cost,len(path)
		if to_use==0:
			r="Search ... " + str(path)
			params['mins']+=1
			return results,True,r
		return results,False,""
	candidates=get_candidates(last_node,fl)
	if len(candidates)==0:
		print "Something really bad"
		sys.exit(1)
	if len(candidates)==1:
		candidate=candidates[0]
		if last_node not in fl:
			fl[last_node]={}
		if candidate not in fl[last_node]:
			fl[last_node][candidate]=0
		fl[last_node][candidate]+=1
		path.append(candidate)
		results.append((len(path),fl,path,to_use+edges[last_node][candidate][0]))
	else:
		for candidate in candidates:
			flr=deepcopy(fl)
			#some quick sanity checks
			if last_node not in flr:
				flr[last_node]={}
			if candidate not in flr[last_node]:
				flr[last_node][candidate]=0
			flr[last_node][candidate]+=1
				
			#lets find out if taking it out gives the same cost
			max_cost,g=print_graph(flr,extra=candidate)
			#print path,candidate
			local_fl,flow_cost=get_flow(g)
			if max_cost==flow_cost:
				if max_cost>max_so_far[0]+100 or max_cost>-100:
					print max_cost,flow_cost,len(path),to_use
					sys.stdout.flush()
					max_so_far[0]=max_cost
				local_path=deepcopy(path)
				local_path.append(candidate)
				results.append((len(local_path),flr,local_path,to_use+edges[last_node][candidate][0]))
			else:
				print >> sys.stderr, "ERROR !!!" , max_cost
				sys.exit(1)
	return results,False,""
	

params['walks']=read_graph(sys.stdin)
print params
fname='/dev/shm/out_%d' % pid
sx,g=print_graph()
flx,cx=get_flow(g)
flx_keys=flx.keys()
random.shuffle(flx_keys)



#try to find the components
k={} #which edge belongs to what component
rk={} #reverse k
cc=0
sz={} #flow per component
for f in flx_keys:
	if f not in k:
		cc+=1
		rk[cc]=set(f)
		sz[cc]=0
		q=[f]
		while len(q)!=0:
			x=q.pop(0)
			new_children=[]
			if x in flx:
				for c in flx[x]:
					if c not in k and flx[x][c]>0:
						k[c]=cc
						rk[cc].add(c)
						new_children.append(c)
						sz[cc]+=flx[x][c]
			q+=new_children

print "There are " ,cc
print sz
#sys.exit(1)	
			


def look(start,end,cc):
	#print_graph(fname)
	#flx,cost=get_flow(fname)
	
	q=[(1,{},[start_node],sz[cc])]
	while len(q)>0:
		tpl=q.pop()
		results,d,r = search(tpl)
		if d or time.time()>end:
			return r
		for result in results:
			if result[2][0]==result[2][-1]:
				print "loop found"
				#now we find a new spot to start from in this cc
				for 
			q.append(result)
	return ""


origin_time = time.time()
max_time=60*60*60

d=60*60*12
while True:
	b=[]
	for cc in sz:
		print "working on cc ",cc
		print >> sys.stderr, "working on cc ",cc
		start = time.time()
		r=look(start,start+d,cc)
		#if time.time()<start+d:
		if len(r)>0:
			#found a solution
			print "FOUND A SOLUTION"
			print >> sys.stderr, "FOUND A SOLUTION"
			b.append(r)
		else:
			#didnt find one
			#d=int(d*1.2)
			print "NO SOLUTION FOUND"
			print >> sys.stderr, "NO SOLUTION FOUND"
		print >> sys.stderr, "Restarting with ", d
	if len(b)>0:
		print "\n".join(b)
		sys.stdout.flush()
