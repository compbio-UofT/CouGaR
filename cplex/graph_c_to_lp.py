

import sys
import gzip


if len(sys.argv)!=8:
	print "%s problemxfile.gz edge_lookup original_edges_wmapq Q somaticQ m output_problem_file" % sys.argv[0]
	sys.exit(1)

pb_fname=sys.argv[1]
el_fname=sys.argv[2]
oe_fname=sys.argv[3]
Q=int(sys.argv[4])
somaticQ=int(sys.argv[5])
m=int(sys.argv[6])
oprob=sys.argv[7]
oprobf=open(oprob,'w')
#Q=500

nodes={}
ins={}
outs={}
edges=[]

'''
Somatic edges types
	0 + + >->
	1 - - <-<
	2 + - >-<
	3 - + <->
'''

edge_lookup=[]
somatic_edges={}
oe={}

def read_oe(fn):
	#chr1:1613433 31.6204
	#chr2:136003458 44.4231
	#chr2:136003458 44.4231
	#chr1:1613433 31.6204
	h=open(fn)
	for line in h:
		line=line.strip().split()
		#print >> sys.stderr, "Enabled chr25 / chrM!"
		#if line[0][:5]=="chr25":
		#	continue
		f=line[0]
		mq=float(line[1])
		#mq=min(35,float(line[1])) #this is for MAPQ 
		oe[f]=mq
	h.close()
	
read_oe(oe_fname)

def get_q(f1,f2):
	sf1=0
	if f1 in oe:
		sf1=oe[f1]
	else:
		print >> sys.stderr, "MISSING ",f1
	sf2=0
	if f2 in oe:
		sf2=oe[f2]	
	else:
		print >> sys.stderr, "MISSING ",f2
	#return m*(1+(35-sf1)/35+(35-sf2)/35)*somaticQ
	return m*(1+(1-sf1)+(1-sf2))*somaticQ


def add_edge(fn,tn,ty,cost,cap):
	eid=len(edges)
	edge_lookup.append((eid,fn,tn,ty))
	edges.append((cost,cap))
	if ty==0:
		outs[fn].append(eid)
		ins[tn].append(eid)
	elif ty==1:
		ins[fn].append(eid)
		outs[tn].append(eid)
	elif ty==2:
		outs[fn].append(eid)
		outs[tn].append(eid)
	elif ty==3:
		ins[fn].append(eid)
		ins[tn].append(eid)


def readxpb(fname):
	h=gzip.open(fname)
	print_arc=0
	arc_cost=0
	for line in h:
		if line[0]=='c':
			if line[2]=='G':
				#c Genomic       chr1:712095     chr1:714492     66      4       2397    736     1844
				#ss << "c Genomic\t" << e.posa.str() << "\t" << e.posb.str()  << "\t" << cost << "\t" << cap << "\t" << e.length() << "\t" << ei.normal << "\t" << ei.tumor << endl;
				line=line.split()
				cost=int(line[4])
				cap=int(line[5])
				l=int(line[6])
				normal=int(line[7])
				tumor=int(line[8])
				fn=nodes[line[2]]
				tn=nodes[line[3]]
				add_edge(fn,tn,0,cost,cap)
				print_arc=2
				arc_cost=cost
				print >> oprobf, line[0]+" "+"\t".join(line[1:])
			elif line[2]=='S':
				#c Somatic       chr1:814892     chrY:11934477   1       23      399     237     33
				#ss << "c Somatic\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << ei.type << "\t" << cost << "\t" << cap << "\t" << ei.normal << "\t" << ei.tumor << endl;
				line=line.split()
				ty=int(line[4])
				cost=int(line[5])
				#cost=int(Q)/2
				#cost=somaticQ
				cost=int(get_q(line[2],line[3]))
				cap=int(line[6])
				normal=int(line[7])
				tumor=int(line[8])
				fn=nodes[line[2]]
				tn=nodes[line[3]]
				if (fn,tn) not in somatic_edges:
					somatic_edges[(fn,tn)]=set()
				somatic_edges[(fn,tn)].add(ty)
				#if line[2]=="chrX:154507287" and line[3]=="chrX:154510610":
				#	print line,fn,tn
				#	sys.exit(1)
				#if fn==257490 and tn==257493:
				#	print line
				#	print len(edges)
				#	sys.exit(1)
				add_edge(fn,tn,ty,cost,cap)
				print_arc=2
				arc_cost=cost
				line[5]=str(cost)
				print >> oprobf, line[0]+" "+"\t".join(line[1:])
			elif line[2]=='N':
				#its a node
				#c NODE 19 chr1:808212
				line=line.split()
				nid=int(line[2])
				if nid%2==0:
					nodes[line[3]]=nid
					ins[nid]=[]
					outs[nid]=[]
				print >> oprobf, " ".join(line)
				print >> oprobf, "n\t"+line[2]+"\t0"
			else:
				print >> oprobf, line,
		elif line[0]=='p':
			print >> oprobf, line,
		elif line[0]=='a':
			#a       1       3       0       0       16
			#a       4       2       0       0       16
			#a       3       7       0       0       0
			#a       8       4       0       0       0
			line=line.split()
			start_node=(int(line[1])+1)/2
			end_node=(int(line[2])+1)/2
			if (start_node==1 and end_node==2) or (start_node==2 and end_node==1):
				line[5]=str(m*Q)
				print >> oprobf, "\t".join(line)
			elif start_node<=2 or end_node<=2:
				line[5]="5"
				print >> oprobf, "\t".join(line)
			elif print_arc>0:
				line[5]=str(arc_cost)
				print >> oprobf, "\t".join(line)
				print_arc-=1
		elif line[0]=='n':
			pass
		
	#add in the funky edges	
	outs[2]=[]
	ins[2]=[]
	ins[4]=[]
	outs[4]=[]
	ins[6]=[]
	outs[6]=[]
	add_edge(2,4,0,m*Q,1000)
	for node in nodes:
		nid=nodes[node]
		if nid>6:
			add_edge(4,nid,0,5,1000)
			add_edge(4,nid,2,5,1000)
			add_edge(nid,2,0,5,1000)
			add_edge(nid,2,3,5,1000)
	h.close()




readxpb(pb_fname)

def outputxlpxsolve():
	objxfunc=[]
	for x in range(len(edges[e])):
			if edges[e][x][0]!=0:
				objxfunc.append("%d %s" % (edges[x][0],eid))
	print "min: "," + ".join(objxfunc),";"
	for n in nodes:
		c=[]
		for x in ins[nodes[n]]:
			c.append(" + f%d" % (x))
		for x in outs[nodes[n]]:
			c.append(" - f%d" % (x)) 
		if len(c)!=0:
			print "".join(c),"= 0",";"
	for e in edges:
		for x in range(len(edges[e])):
			eid="%sx%d" % (e,x)
			print " %s <= %d" % (eid,edges[eid][x][1]),';'
	print "int "
	o=[]
	for e in edges:
		for x in range(len(edges[e])):
			eid="%sx%d" % (e,x)
			o.append(eid)
	print ",".join(o),';'

def outputxcplex():
	#now output the lp file
	objxfunc=[]
	for x in range(len(edges)):
		if edges[x][0]!=0:
			objxfunc.append("%d f%d" % (edges[x][0],x))
	print "Minimize\n obj:",
	ll=5
	for x in objxfunc:
		if ll+len(x)+3>=80:
			ll=0
			print ""
			print " +",x,
			ll+=len(x)+3
		else:
			print " +",x,
			ll+=len(x)+3
	print ""
	print "Subject To"
	for n in nodes:
		c=[]
		for x in ins[nodes[n]]:
			c.append(" + f%d" % (x))
		for x in outs[nodes[n]]:
			c.append(" - f%d" % (x)) 
		if len(c)!=0:
			print " n%d: " % nodes[n],"".join(c)," = 0"
			#print "".join(c),"= 0",";"
	#put in the source and sink
	
	print "Bounds"
	for x in range(len(edges)):
		print "0 <= f%d <= %d" % (x,edges[x][1])
	#print "Integers "
	#o=[]
	#for x in range(len(edges)):
	#	o.append("f%d" % x)
	#print " ".join(o)
	print "End"

outputxcplex()

h=gzip.open(el_fname,'w')
h.write("\n".join(map(lambda x : str(x) , edge_lookup)))
h.close()
