


cancers={}
patient={}
lookup={}

cols=4

def get_set(l):
	s=set()
	l=l.replace('-','').split(',')
	for g in l:
		if len(g)>0:
			s.add(g)
	return s

import sys
#TCGA-ER-A197    SKCM    gfloopamps      -       -       -       -
#populate per patient per column
for line in sys.stdin:
	line=line.strip().split()
	if line[0] not in patient:
		patient[line[0]]=[]
		for x in range(cols):
			patient[line[0]].append(set())
	for x in range(cols):
		patient[line[0]][x].update(get_set(line[3+x]))
	# add to lookup
	if line[1] not in cancers:
		cancers[line[1]]=set()
	cancers[line[1]].add(line[0])		
	line[1]='all'
	if line[1] not in cancers:
		cancers[line[1]]=set()
	cancers[line[1]].add(line[0])		
	
	 
#get stats per column per cancer
gl=[]
for x in range(cols):
	gl.append({})
for c in cancers:
	print c+"\t",
	gs=[]
	for x in range(cols):
		gs.append({})
	#build dictionary
	for p in cancers[c]:
		for x in range(cols):
			for g in patient[p][x]:
				if g not in gs[x]:
					gs[x][g]=0
				if g not in gl[x]:
					gl[x][g]=0
				gs[x][g]+=1
				gl[x][g]+=1
	o=[]
	for x in range(cols):
		l=[]
		oo=[]
		ooo=[]
		for g in gs[x]:
			l.append((gs[x][g],g))
		l.sort(reverse=True)
		if len(l)>0:
			mx=l[0][0]
			if mx<=2:
				oo.append('-')
				o.append(",".join(oo))
			else:
				for c,g in l:
					if c==mx:
						oo.append(g)
					elif c==mx-1:
						ooo.append(g)
				o.append(str(mx)+"|"+",".join(oo))
				#if len(ooo)>0:
				#	o[-1]+=("/" +str(mx-1)+"|"+",".join(ooo))
		else:
			oo.append('-')
			o.append(",".join(oo))
	print "\t".join(o)	
