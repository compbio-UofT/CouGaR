#!/usr/bin/python


import sys
from copy import deepcopy

#WARNING MAKING ASSUMPTIONS ABOUT EDGES.......
#BAD ASSUMPTIONS

node_to_pos={}
genomic_edges={}
somatic_edges={}
all_edges={}

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
	return int(s)

if len(sys.argv)!=5:
	print "%s problem_file paths scale output_folder" % sys.argv[0]
	sys.exit(1)



def read_problem_file(filename):
	h=open(filename,'r')
	for line in h:
		line=line.split()
		if line[1]=='NODE':
			#node line
			#c NODE 24 chr1:202709642
			n=int(line[2])
			node_to_pos[n]=to_pos(line[3])
		elif line[1]=='Somatic':
			#Somatic line
			#c Somatic       chr1:202869943  chr1:204056467  2
			f=to_pos(line[2])
			t=to_pos(line[3])
			if f>t:
				print >> sys.stderr, "ERROR" 
				sys.exit(1)
			typ=int(line[4])
			cost=int(line[5])
			if (f,t) not in somatic_edges:
				somatic_edges[(f,t)]=[]
			somatic_edges[(f,t)].append((cost,False,typ))
			somatic_edges[(f,t)].append((cost,False,typ))
		elif line[1]=='Genomic':
			#genomid line
			#c Genomic       chr2:8917257    chr2:8921728
			f=to_pos(line[2])
			t=to_pos(line[3])
			if f>t:
				print >> sys.stderr, "ERROR" 
				sys.exit(1)
			cost=int(line[4])
			if (f,t) not in genomic_edges:
				genomic_edges[(f,t)]=[]
			genomic_edges[(f,t)].append((cost,True,0))
			genomic_edges[(f,t)].append((cost,True,0))
	h.close()
	for (f,t) in somatic_edges:
		if (f,t) not in all_edges:
			all_edges[(f,t)]=[]
		all_edges[(f,t)]+=somatic_edges[(f,t)]
		all_edges[(f,t)].sort()
	for (f,t) in genomic_edges:
		if (f,t) not in all_edges:
			all_edges[(f,t)]=[]
		all_edges[(f,t)]+=genomic_edges[(f,t)]
		all_edges[(f,t)].sort()
	
def annotate_path(p):	
	used={}
	path=[]
	for x in range(1,len(p)):
		f=node_to_pos[p[x-1]]
		t=node_to_pos[p[x]]
		# not the source and sink chromosomes
		if f[0]>0 and t[0]>0:
			#get the canonical
			c=(f,t)
			if f>t:
				c=(t,f)
			if c not in all_edges:
				print >> sys.stderr, " Dropping edge " ,f,t
				sys.exit(1)
			else:
				if c not in used:
					used[c]={'genomic':0,'somatic':0}
				#find out if we are using somatic or genomic
				genomic=all_edges[c][used[c]['genomic']+used[c]['somatic']][1]	
				#print len(all_edges[c]),used[c],c
				if genomic:
					path.append((f,t,True))
					used[c]['genomic']+=1
				else:
					path.append((f,t,False))
					used[c]['somatic']+=1
	return used,path



		
def pr(pf,pt,pl,paths):
	d={}
	s=0
	for path in paths:
		last=-1
		added=set()
		#does the forward looking pass
		for x in range(1,len(path)):
			f=path[x-1]
			t=path[x]
			if f<=4 or t<=4:
				#just the top start nodes
				continue
			elif f<=6 or t<=6:
				#hit an end node
				break
			if pf==f and pt==t:
				#found one
				last=0
			else:
				if last>=0:
					last+=1
					if last==pl:
						z=tuple(path[x-1-last:x+1])
						if z not in d:
							d[z]=0
						if z not in added:
							d[z]+=1
							added.add(z)
						#d[z]+=1
						#s+=1
		s+=1
	print s,d
	dist={}
	mx=0
	for z in d:
		dist[z]=d[z]/float(s)
		mx=max(mx,dist[z])
	return dist,mx

def reduce_solid_paths(sp):
	to_remove=set()
	for p in sp:
		#check start with genomic
		f=node_to_pos[p[0]]
		t=node_to_pos[p[1]]
		if (f,t) not in genomic_edges and (t,f) not in genomic_edges:
			print >> sys.stderr, "X Dropping path because does not start with genomic edge ", p
			to_remove.add(p)	
			continue
		#check end with genomic
		f=node_to_pos[p[-2]]
		t=node_to_pos[p[-1]]
		if (f,t) not in genomic_edges and (t,f) not in genomic_edges:
			print >> sys.stderr, "X Dropping path because does not end with genomic edge ", p
			to_remove.add(p)	
			continue
		for x in range(1,3):
			bwd=tuple(p[:-x])
			fwd=tuple(p[x:])
			if fwd in sp:
				to_remove.add(fwd)
			if bwd in sp:
				to_remove.add(bwd)
	r=set()
	for p in sp:
		if p not in to_remove:
			r.add(p)
	print >> sys.stderr, "Original size: %d, now %d" % (len(sp),len(r))
	return r
	

def collapse_path(p):
	path=[]
	for f,t,g in p:
		if len(path)>0:
			if g and path[-1][2]:
				#this should be def mean that we can collapse
				if path[-1][1]!=f:
					print >> sys.stderr, " MAJOR FAIL"
					sys.exit(1)
				path[-1][1]=t 
			else:
				path.append([f,t,g])
		else:
			path.append([f,t,g])
	return path
		

def pos_to_str(p):
	x=""
	if p[0]==23:
		x="X"
	elif p[0]==24:
		x="Y"
	else:
		x=str(p[0])
	return "chr"+x,str(p[1])
					
def e_to_bed(f,t):
	if f>t:
		z=f
		f=t
		t=z
	chra,coorda=pos_to_str(f)
	chrb,coordb=pos_to_str(t)
	return chra+"\t"+coorda+"\t"+coordb
	

def cpath_to_bed(cp):
	r=[]
	for f,t,g in cp:
		if g:
			r.append(e_to_bed(f,t))
	return "\n".join(r)+"\n"
		
		
	
def read_paths_file(filename,output_folder):
	paths=set()
	edges=set()
	h=open(filename,'r')
	for line in h:
		if line.find('Search ... ')==0:
			p=None
			try:
				p=tuple(eval(line.replace('Search ... ','')))
			except:
				pass
			if p!=None:
				paths.add(p)
				e=set()
				for x in range(1,len(p)):
					f=p[x-1]
					t=p[x]
					if f<=6 or t<=6:
						continue
					e.add((f,t))
				#edges=edges.union(e)
				edges.update(e)
				if len(paths)>1000:
					print >> sys.stderr, " ONLYF IRST 50 paths"
					break
	h.close()


	solid_paths=set()

	#base settings
	slack=1
	#edges=set([(116, 94)])
	#edges=set([(63,65)])
	#edges=set([(98,96)])
	#now for each edge
	for f,t in edges:
		l=1
		p_cutoff=0.999
		x=slack
		#last=set()
		while x>0:
			dist,mx=pr(f,t,l,paths)
			#print dist
			if mx<p_cutoff:
				x-=1
			else:
				x=slack
				for z in dist:
					if dist[z]>=p_cutoff:
						solid_paths.add(z)
				#last=dist.keys()
			l+=1
		#print last
		#print last
		#solid_paths=solid_paths.union(last)
		#solid_paths.update(last)
	solid_paths=reduce_solid_paths(solid_paths)
	path_id=0
	for path in solid_paths:
		path_id+=1
		used,path=annotate_path(path)
		cpath=collapse_path(path)
		if len(cpath)>1:
			#print "# INFO" + filename + " solid_path_len" , len(path)
			#print "# NODEPATH" , path
			#print "# PATH " ,  path
			#print "# COLLAPSEDPATH " ,  cpath
			#write collapsed path to file
			h=open(output_folder+"/"+str(path_id)+"_cpath.bed",'w')
			h.write(cpath_to_bed(cpath))
			h.close()
			#write out the graph file
			sout=""
			l_genomic_edges=[]
			for f,t in genomic_edges:
				u=0
				if (f,t) in used:
					u=used[(f,t)]['genomic']*scale
				l_genomic_edges.append((f,t,0,u))	
			sout+=str(l_genomic_edges)+"\n"
			l_somatic_edges=[]
			for f,t in somatic_edges:
				u=0
				if (f,t) in used:
					u=used[(f,t)]['somatic']*scale
				l_somatic_edges.append((f,t,somatic_edges[(f,t)][0][2],u))
			sout+=str(l_somatic_edges)+"\n"
			sout.replace(',','')
			#write it out to file
			h=open(output_folder+"/"+str(path_id)+"_graph",'w')
			h.write(sout)
			h.close()

problem_filename=sys.argv[1]
paths_filename=sys.argv[2]
scale=float(sys.argv[3])
output_folder=sys.argv[4]

read_problem_file(problem_filename)
read_paths_file(paths_filename,output_folder)

