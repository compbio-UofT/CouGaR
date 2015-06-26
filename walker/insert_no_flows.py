

import sys


if len(sys.argv)!=3:
	print sys.argv[0],"g_lp hmm"
	sys.exit(1)

g_lp_fn=sys.argv[1]
hmm_fn=sys.argv[2]

somatic=[]
genomic=[]


#read in the g_lp
h=open(g_lp_fn)
line=h.readline()
title="#"
print >> sys.stderr, "READING" , g_lp_fn
while line[0]=='#':
	title=line.strip()
	line=h.readline()
genomic=eval(line.strip().replace(' ',','))
line=h.readline()
somatic=eval(line.strip().replace(' ',','))

def score(l):
	su=0
	st=0
	for (s,e,x,u,n,t) in l:
		su+=abs(t-u*n)
		st+=abs(t-2*n)
	return su,st
h.close()

def chr_to_n(c):
	c=c.lower()
	if c[:3]=='chr':
		c=c[3:]
	if c[0]=='x':
		return 23
	if c[0]=='y':
		return 24
	if c[0]=='m':
		#print >> sys.stderr, "ENABLED CHR25!"
		return 25
	return int(c)

#read in the hmm
#						T	N
#2      chr12:25606866 chr12:25609266    2400    632     418
#1      chr12:25609266 chr12:25611666    2400    506     475
#2      chr12:25611666 chr12:25657266    45600   14435   7270                  
#1      chr12:25657266 chr12:25659666    2400    574     441
#2      chr12:25659666 chr12:25662066    2400    893     445
hmm_edges=[]
min_cp=5
min_len=40000
h=open(hmm_fn)
for line in h:
	if line[0]=='#':
		continue
	line=line.strip().split()
	cp=int(line[0])
	if cp<min_cp:
		continue	
	s=(chr_to_n(line[1].split(":")[0]),int(line[1].split(":")[1]))
	e=(chr_to_n(line[2].split(":")[0]),int(line[2].split(":")[1]))
	
	#ok good length and good copy count
	if len(hmm_edges)>1:
		pdoc=float(hmm_edges[-1][5])/(hmm_edges[-1][4]+1)
		cdoc=float(line[4])/(int(line[5])+1)
		if hmm_edges[-1][1][0]==s[0] and ( (e[1]-s[1]<3000 or abs(pdoc-cdoc)<3) and hmm_edges[-1][1]==s):
			#print pdoc,cdoc,s,e
			hmm_edges[-1]=(hmm_edges[-1][0],e,0,0,hmm_edges[-1][4]+int(line[5]),hmm_edges[-1][5]+int(line[4]))
		else:
			hmm_edges.append((s , e, 0, 0, int(line[5]), int(line[4])))
	else:
		hmm_edges.append((s , e, 0, 0, int(line[5]), int(line[4])))
h.close()

new_edges=[]
for edge in hmm_edges:
	if edge[1][1]-edge[0][1]>=min_len:
		new_edges.append(edge)

hmm_edges=new_edges


x=[((0,0),(0,0),0,0,0,0)] + genomic[:] + [((25,0),(25,0),0,0,0,0)]

new_edges=[]
for edge in hmm_edges:
	for i in range(1,len(x)):
		c=x[i]
		p=x[i-1]
		#totally diff chr
		if p[1]==c[0]:
			continue
		interval=(p[1],c[0])
		#find out overlap with this interval
		s=max(interval[0],edge[0])
		e=min(interval[1],edge[1])
		if e>s:
			new_edges.append((s,e,0,0,edge[4],edge[5]))

		
gsu,gst=score(genomic)
	
genomic+=new_edges
genomic.sort()

gsup,gstp=score(genomic)
print >> sys.stderr, float(gsu)/(gst+0.01), float(gsup)/(gstp+0.01)

print title + " with HMM edges mincp" , min_cp, "minlen", min_len
print str(genomic).replace(',','')
print str(somatic).replace(',','')
		
	
