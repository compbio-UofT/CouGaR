import sys

from random import choice

if len(sys.argv)!=4:
	print "%s flow_file out out_master"
	sys.exit(1)

flow_filename=sys.argv[1]
out_filename=sys.argv[2]
out_master_filename=sys.argv[3]



flow_file=open(flow_filename)

edges={}
ins={}
outs={}

nodes=set()
funky_nodes=set()

for line in flow_file:
	line=line.strip().split()
	if line[0]=='f':
		fl=int(line[3])
		if fl==0:
			continue
		fn=int(line[1])
		tn=int(line[2])
		if fn not in edges:
			edges[fn]={}
		if tn not in edges[fn]:
			edges[fn][tn]=0
		edges[fn][tn]+=fl
		if tn not in ins:
			ins[tn]=0
		ins[tn]+=fl
		if fn not in outs:
			outs[fn]=0
		outs[fn]+=fl
		nodes.add(fn)
		nodes.add(tn)
		if abs(fn-tn)!=2:
			funky_nodes.add(fn)

for x in ins:
	if ins[x]!=outs[x]:
		print >> sys.stderr, "FAIL",x
		sys.exit(1)

for x in outs:
	if ins[x]!=outs[x]:
		print >> sys.stderr, "FAIL",x
		sys.exit(1)
		

def get_start_edge(e):
	es=[]
	for fn in e:
		for tn in e[fn]:
			if e[fn][tn]<=0:
				print >> sys.stderr, "FAILDE!"
				sys.exit(1)
			for x in range(e[fn][tn]):
				es.append((fn,tn))
	if len(es)==0:
		return False
	return choice(es)			

def decrease(e,fn,tn):
	e[fn][tn]-=1
	if e[fn][tn]==0:
		e[fn].pop(tn)
	if len(e[fn])==0:
		e.pop(fn)



def reverse_path(p):
	p=list(p)
	for x in range(len(p)):
		if p[x]%2==0:
			p[x]-=1
		else:
			p[x]+=1
	p.reverse()
	return tuple(p)
	

def canonical_loop(l):
	mn=min(l)
	if mn%2==0:
		l=reverse_path(l)
	l=l[l.index(min(l)):]+l[:l.index(min(l))] #bring min to front
	return tuple(l)


#given loop and subseq s, give a histogram of things to come next
def ocount(l,s):
	c={}
	sz=len(l)
	szz=len(s)
	for i in range(sz):
		missed=False
		for x in range(szz):
			if l[(i+x)%sz]!=s[x]:
				missed=True
				break
		if not missed:
			z=l[(i+szz)%sz]
			if z!=None:
				if z not in c:
					c[z]=0
				c[z]+=1	
	return c

def is_subsequence(s,m):
	sz=len(m)
	szz=len(s)
	for i in range(sz-szz+1):
		missed=False
		for x in range(szz):
			if m[i+x]!=s[x]:
				missed=True
				break
		if not missed:
			return True
	return False


def print_seq(ls):
	if len(ls)==0:
		return ""
	os=[]
	for s in ls:	
		if len(s)==1:
			os.append(str[0])
			continue
		o=""
		d=None
		l=0
		for x in range(len(s)-1):
			if s[x]-s[x+1]!=d:
				if len(o)>0:
					if l>0:
						o+="-"+str(s[x])
					else:
						o+=","+str(s[x])
					#o+=","+str(s[x+1])
				else:
					o=str(s[x])
				d=s[x]-s[x+1]
				l=0
			else:
				l+=1
		o+="-"+str(s[-1])
		os.append(o)
	return "|".join(os)
		

def is_primary(l):
	l=canonical_loop(l)
	try:
		sa=l.index(1)
		return True
	except:
		pass
	return False	

def loop_to_loops(l):
	l=canonical_loop(l)
	sa=-1
	sb=-1
	#l cannot have repeat element, but can go through source/sink on other side
	try:
		sa=l.index(1)
	except ValueError:
		pass
	try:
		sb=l.index(4)
	except ValueError:
		pass
	if sa>=0 and sb>=0:
		#this is a double loop
		l1=l[:sb]
		l2=canonical_loop(l[sb:])
		return [(tuple(l1),),(tuple(l2),)]
	else:
		return [(tuple(l),)]

all_loops={}


def join_loops(l1,l2):
	nl=set()
	for l in l1:
		nl.add(l)
	for l in l2:
		nl.add(l)
	nl=list(nl)
	nl.sort()
	return tuple(nl)

def loop_to_loop(l1,l2):
	if len(l2)!=1:
		print >> sys.stderr, "FAILEDS AT THIS"
		sys.exit(1)
	l2=l2[0]
	#check if the loop is already in l1
	for loop in l1:
		if loop==l2:
			return ()
	#check for overlaping node
	for loop in l1:
		for x in loop:
			if x in l2:
				#they overlap
				l3=list(l1)+[l2]
				l3.sort()
				return tuple(l3)
	return ()
	
		

#take in contigs
# { '(1,3,4,16,347,17,547,245,724572)':count ...}
#same for loops
def cplexout(lines,loops,candidates,filename):
	#should really remove length 1 and copy count 1 candidates
	#lmbda=1
	#lmbda=-0.3
	lmbda=60
	new_candidates=[]
	for c in candidates:
		r=[]
		q=[]
		for x in c:
			r.append(x)
			q.append(c[x])
		r=tuple(r)
		q=tuple(q)
		new_candidates.append((r,q))
	candidates=new_candidates
	
	of=open(filename,'w')
	#index the required and candidates
	idxs={}
	line_idxs=[]
	idx_to_line={}
	loop_idxs=[]
	idx_to_loop={}
	candidate_idxs=[]
	idx_to_candidate={}
	for x in lines:
		line_idxs.append(len(idxs))
		idx_to_line[len(idxs)]=x
		idxs[x]=len(idxs)
	for x in loops:
		loop_idxs.append(len(idxs))
		idx_to_loop[len(idxs)]=x
		idxs[x]=len(idxs)
	for x in candidates:
		candidate_idxs.append(len(idxs))
		idx_to_candidate[len(idxs)]=x
		idxs[x]=len(idxs)
	ks=idxs.keys()
	#generate the objective 
	obj=[]
	for x in range(len(idxs)):
		obj.append(str(lmbda)+' m'+str(x))
		obj.append('i'+str(x))
	#for x in idxs:
	#	obj.append(str(lmbda)+' m'+str(idxs[x]))
	#for x in idxs:
	#	obj.append('i'+str(idxs[x]))
	print >> of, "Minimize\n obj:"," + ".join(obj)
	
	#print the conditions
	print >> of, "Subject To"
	for x in line_idxs:
		l=idx_to_line[x]
		terms=[]
		terms.append("m"+str(x))
		for y in candidate_idxs:
			c,cp=idx_to_candidate[y]
			if l[0] in c:
				# print M * (how many times c has l in it)
				i=c.index(l[0])
				terms.append(str(cp[i]) + " m"+str(y))
		print >> of,  " l"+str(idxs[l])+":\t" + " + ".join(terms) + " = " + str(lines[l])
	for x in loop_idxs:
		l=idx_to_loop[x]
		terms=[]
		terms.append("m"+str(x))
		for y in candidate_idxs:
			c,cp=idx_to_candidate[y]
			if l[0] in c:
				# print M * (how many times c has l in it)
				i=c.index(l[0])
				terms.append(str(cp[i]) + " m"+str(y))
		print >> of, " o"+str(idxs[l])+":\t" + " + ".join(terms) + " = " + str(loops[l])
	for x in range(len(idxs)):
		if x in idx_to_line:
			print >> of , "\\ LINE " , print_seq(idx_to_line[x])	
		elif x in idx_to_loop:
			print >> of , "\\ LOOP " , print_seq(idx_to_loop[x])	
		elif x in idx_to_candidate:
			print >> of , "\\ CANDIDATE " , print_seq(idx_to_candidate[x][0]), idx_to_candidate[x][1]	
		else:
			print "MAJOR ERROR Xsadf"
			sys.exit(1)
		print >> of,  " lx"+str(x)+":\t" + " m"+ str(x) +" - 1000 i"+str(x)+ " < 0"
		

	#print the bounds
	print >> of, "Bounds"
	for b in range(len(idxs)):
		print >> of, "0 <= m" + str(b) + " <= 1000"

	print >> of, "Binary"
	print >> of, "\t"+ " ".join(map(lambda x : "i"+str(x),range(len(idxs))))
		
	print >> of, "General"
	print >> of, "\t"+ " ".join(map(lambda x : "m"+str(x),range(len(idxs))))

	print >> of , "End"
	of.close()



def common_factor(depth,s):
	for x in range(2,depth):
		good=True
		for loop in s:
			if s[loop]%x!=0:
				good=False
				break
		if good:
			return True
	return False
	

def sloop_candidates(depth,s):
	candidates=sloop_candidates_r(depth,s)
	new_candidates=[]
	for x in candidates:
		if common_factor(depth,x):
			continue
		new_candidates.append(x)
	return candidates

def sloop_candidates_r(depth,s):
	r=[]
	if len(s)==1:
		for x in range(depth):
			r.append({s[0]:(x+1)})
	else:
		qq=sloop_candidates(depth,s[1:])
		for q in qq:
			for x in range(depth):
				d=q.copy()
				d[s[0]]=x+1
				r.append(d)
	return r


def decompose(oe,lines,loops):
	e=edges_copy(oe)
	sz=0
	for fn in e:
		for tn in e[fn]:
			sz+=e[fn][tn]	
	
	#find the neighbours
	stack=[]
	while True:
		if len(stack)<2:
			start_edge=get_start_edge(e)
			if start_edge==False:
				break
			decrease(e,start_edge[0],start_edge[1])
			sz-=1
			stack=[start_edge[0],start_edge[1]]
		fn=stack[-1]
		if len(e[fn])==1:
			tn=e[fn].keys()[0]
			decrease(e,fn,tn)
		else:
			ns=[]
			for tn in e[fn]:
				for x in range(e[fn][tn]):
					ns.append(tn)
			if len(ns)==0:
				break
			#choose a neighbour 
			tn=choice(ns)
			decrease(e,fn,tn)
		sz-=1
		if tn in stack:
			#remove the loop!
			idx=stack.index(tn)
			#loop=canonical_loop(stack[idx:])
			new_loops=loop_to_loops(stack[idx:])
			for loop in new_loops:
				print print_seq(loop),
				d=None
				if is_primary(loop[0]):
					print "LINE"
					d=lines
				else:
					print "LOOP"
					d=loops
				if loop not in d:
					d[loop]=0
				d[loop]+=1
			stack=stack[:idx]
			#if len(loops)>3 and len(lines)>3:
			#	break
		stack.append(tn)


	return 

def reachable_without(oe,fn,tn):
	e=edges_copy(oe)
	decrease(e,fn,tn)
	visited=set([1,2])
	if tn in (1,2):
		return True
	q=[1,2]
	while len(q)>0:
		n=q.pop() #DFS if pop from back, BFS if pop from front?
		for t in e[n]:
			if tn==t:
				return True
			if t not in visited:
				visited.add(t)
				q.append(t)
	return False
	
		
	

def connected_components_without(oe,fn,tn):
	e=edges_copy(oe)
	decrease(e,fn,tn)
	return unsafe_connected_components(e)

def unsafe_connected_components(e):
	cc=0
	while len(e)>0:
		#grab the first node
		stack=[e.keys()[0]]
		while len(stack)>0:
			fn=stack.pop()
			if fn in e:
				for tn in e[fn]:
					stack.append(tn)
				e.pop(fn)
		cc+=1
	return cc

def connected_components(oe):
	e=edges_copy(oe)
	return unsafe_connected_components(e)

def edges_copy(e):
	ne={}
	for x in e:
		ne[x]=e[x].copy()
	return ne

out_file=open(out_filename,'w')
out_master_file=open(out_master_filename,'w')
lines={}
loops={}
for x in range(100):
	decompose(edges,lines,loops)

for x in lines:
	print >> out_file, str(x)
	print >> out_master_file, "LINE", str(x)

for x in loops:
	print >> out_file, str(x)
	print >> out_master_file, "LOOP", str(x)



