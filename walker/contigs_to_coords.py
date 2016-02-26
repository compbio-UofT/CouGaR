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
			#print >> sys.stderr, "CHRM is enabeled"
			return 25
	return int(s)

if len(sys.argv)!=9:
	print "%s problem_file contig_file ip_sol m genes_file oncogenes name groups" % sys.argv[0]
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


#read in the integer programming solution
#return a list of tuples, (contigID , copy count)
def read_sol_file(fname):	
	ip_sol_file=open(fname)
	cis=[]
	for line in ip_sol_file:
		if line[0]=="c":
			line=line.split()
			if int(line[1])>0:
				cis.append((int(line[0][1:]),int(line[1])*m))
	return cis

#takes a file name of decomposed paths and returns a list of tuples, (type,path)
def read_simple_contigs_file(fname):
	simple_contigs=[]
	simple_contigs_f=open(fname)
	for line in simple_contigs_f:
		if line[0]!="L": #make sure this is a line or a loop
			continue
		ty=line[:4] #the type, line or loop
		line=line[5:] #the description of the pay
		simple_contigs.append((ty,eval(line.strip())[0])) #append the type and the list of nodes
	return simple_contigs

def read_genes(fname):
	f=open(fname)
	#((1 10003485) (1 10045556) NMNAT1 +)
	#((1 100111430) (1 100160097) PALMD +)
	genes=[]
	for line in f:
		if line[0]=="#":
			continue
		line=line.replace("(",'').replace(")","").strip().split()
		chr=int(line[0])
		s=int(line[1])
		e=int(line[3])
		if e<s:
			print "BIG ERROR"
			sys.exit(1)
		g=line[4]
		st=line[5]
		genes.append((chr,s,e,g,st))
	return genes
	

def read_onco(onco_filename):
	onco=set()
	f=open(onco_filename)
	for line in f:
		if line[0]=="#":
			continue
		onco.add(line.strip())
	return onco

#take a list of genes a genomic interval
#port a dictionary of gene occurences and gene interrupts 
def gene_intersect(ogenes,z):
	chr,s,e=z
	genes={}
	genes_interrupted={}
	for gchr,gs,ge,gg,gst in ogenes: # gene chromosome, gene start, gene end, gene name, gene strand
		if gchr!=chr:
			continue
		mx=max(gs,s)
		mn=min(ge,e)
		if mx<mn:
			#some overlap
			if ge<=e and gs>=s:
				#contained
				#genes.add(gg)
				if gg not in genes:
					genes[gg]=0
				genes[gg]+=1
			else:
				#partial
				#genes_interrupted.add(gg)
				if gg not in genes_interrupted:
					genes_interrupted[gg]=0
				genes_interrupted[gg]+=1
	return genes,genes_interrupted
		

		


problem_filename=sys.argv[1]
simple_contigs_filename=sys.argv[2]
ip_sol_filename=sys.argv[3]
m=int(sys.argv[4])
genes_filename=sys.argv[5]
onco_filename=sys.argv[6]
name=sys.argv[7]
groups_filename=sys.argv[8]

groups={}
groups_f=open(groups_filename)
for line in groups_f:
	line=line.strip().split()
	groups[line[0]]=line[1]


genes=read_genes(genes_filename)
cis=read_sol_file(ip_sol_filename) # cis a list of tuples, (contigID,copyCOunt)
simple_contigs=read_simple_contigs_file(simple_contigs_filename)
read_problem_file(problem_filename)
onco=read_onco(onco_filename)
simple_walks=[]

print cis
#map the simple walks from nodes back to genomic intervals
#simple contigs is tuples of (type, path) - type in ['LOOP','LINE'], path is a list of nodes
#simple walks becomes a list of tuples (type in ['LOOP','LINE'],copycount, list of genomic edges)
for ci,cm in cis: #ci - contig id, cm - contig copy count (contig mupliple)
	ty,c=simple_contigs[ci]
	walk=[]
	for x in range(len(c)):
		fn=c[x] #from node
		tn=c[(x+1)%len(c)] #to node
		fnp=node_to_pos[fn]
		tnp=node_to_pos[tn]
		#if its not in the dicitonary then this is a somatic edge?
		if (fnp,tnp) not in genomic_edges and (tnp,fnp) not in genomic_edges:
			continue
		#if this is not a sequential edge then just add it, otherwise extend
		if len(walk)==0 or walk[-1][1]!=fnp:
			walk.append([fnp,tnp])
		else:
			walk[-1][1]=tnp
	print walk
	simple_walks.append((ty,cm,walk))

segs={}
gf_lines={}
gi_lines={}
gf_loops={}
gi_loops={}
lines=[]
loops=[]
total_new_dna=0
i=0
for ty,cm,w in simple_walks:
	i+=1
	print "WALK",i
	l=0
	lsegs={}
	for fn,tn in w:
		schr,scoord=fn
		echr,ecoord=tn
		if schr!=echr:
			print "NFDOFN"
			sys.exit(1)
		if scoord>ecoord:
			t=scoord
			scoord=ecoord
			ecoord=t
		if scoord>ecoord:
			print "SDFS"
			sys.exit(1)
		
		#update segs
		k=(schr,scoord,ecoord)
		if k not in lsegs:
			lsegs[k]=0
		lsegs[k]+=1

		#update length
		l+=ecoord-scoord

		#update genes and interrupted genes		
		g1,g2=gene_intersect(genes,(schr,scoord,ecoord))
		if ty=="LOOP":
			for g in g1:
				if g not in gf_loops:
					gf_loops[g]=0
				gf_loops[g]+=cm
			for g in g2:
				if g not in gi_loops:
					gi_loops[g]=0
				gi_loops[g]+=cm
		elif ty=="LINE":
			for g in g1:
				if g not in gf_lines:
					gf_lines[g]=0
				gf_lines[g]+=cm
			for g in g2:
				if g not in gi_lines:
					gi_lines[g]=0
				gi_lines[g]+=cm
		else:
			print "BIG ERRPR 2"
			sys.exit(1)
	#updates
	total_new_dna+=l*cm
	if ty=="LOOP":
		loops.append((cm,l))
	elif ty=="LINE":
		lines.append((cm,l))
	else:
		print "BIGBIG ERROR"
		sys.exit(1)
	for schr,scoord,ecoord in lsegs:
		if schr not in segs:
			segs[schr]={}
		k=cm*lsegs[(schr,scoord,ecoord)]
		if k not in segs[schr]:
			segs[schr][k]=0
		segs[schr][k]+=1
	print ty,cm,l,len(w),w


def bin_genes(bins,genes):
	genes_amp={}
	for bin in bins:
		genes_amp[bin]={}
	for gene in genes:
		for b in bins:
			if genes[gene]>=b:
				genes_amp[b][gene]=genes[gene]
				#break
	#format string
	o=[]
	for bin in bins:
		oo=[]
		oncos=[]
		for gene in genes_amp[bin]:
			if gene in onco:
				oo.append(gene+"*")
				oncos.append(gene+"*")
			else:
				oo.append(gene)
				pass	
		if len(oo)==0:
			oo.append('-')
		if len(oo)>6000000:
			o.append(",".join(oncos+["..."]))
		else:
			o.append(",".join(oo))
	return "\t".join(o)
	

def bin_contigs(c,sz,bins):
	contigs_amp={}
	for bin in bins:
		contigs_amp[bin]=[0,0,0,0]
	for cm,l in c:
		for b in bins:
			if cm>=b:
				if l>sz:
					contigs_amp[b][2]+=1
					contigs_amp[b][3]+=l
				else:
					contigs_amp[b][0]+=1
					contigs_amp[b][1]+=l
				break
	o=[]
	for bin in bins:
		o.append("\t".join(map(str,contigs_amp[bin])))
	return "\t".join(o)


def bin_segs(segs,bins):
	segs_bin={}
	for bin in bins:
		segs_bin[bin]={}
	for chr in segs:
		for x in segs[chr]:
			for bin in bins:
				if x>=bin:
					if chr not in segs_bin[bin]:
						segs_bin[bin][chr]=0
					segs_bin[bin][chr]+=segs[chr][x]
					#break
	o=[]
	for bin in bins:
		oo=[]
		for chr in segs_bin[bin]:
			if chr==23:
				oo.append("X/%d" % (segs_bin[bin][chr]))
			elif chr==24:
				oo.append("Y/%d" % (segs_bin[bin][chr]))
			elif chr==25:
				oo.append("M/%d" % (segs_bin[bin][chr]))
			else:
				oo.append("%d/%d" % (chr,segs_bin[bin][chr]))
		if len(oo)==0:
			oo.append('-')
		o.append(",".join(oo))
	return "\t".join(o)

#compute the gene stats
amp_bins=[6,15,30,72]
amp_bins.sort(reverse=True)

gf_loop_amps=bin_genes(amp_bins,gf_loops)
gi_loop_amps=bin_genes(amp_bins,gi_loops)
gf_line_amps=bin_genes(amp_bins,gf_lines)
gi_line_amps=bin_genes(amp_bins,gi_lines)

print name + "\t" + groups[name] + "\tgfloopamps\t" + gf_loop_amps
print name + "\t" + groups[name] + "\tgflineamps\t" + gf_line_amps

print name + "\t" + groups[name] + "\tgiloopamps\t" + gi_loop_amps
print name + "\t" + groups[name] + "\tgilineamps\t" + gi_line_amps

#compute contig stats
amp_bins=[0,6,15,30,72]
amp_bins.sort(reverse=True)


print name + "\t" + groups[name] + "\tloopamps\t" + bin_contigs(loops,500000,amp_bins)
print name + "\t" + groups[name] + "\tlineamps\t" + bin_contigs(lines,500000,amp_bins)


#compute the segs
amp_bins=[3,6,15,30,72]
amp_bins.sort(reverse=True)
print name + "\t" + groups[name] + "\tsegs\t" + bin_segs(segs,amp_bins)

