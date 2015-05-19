


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

genes_vs_tumors=[]
tumors_vs_genes=[]
for x in range(cols):
	genes_vs_tumors.append({})
	tumors_vs_genes.append({})

#TCGA-ER-A197    SKCM    gfloopamps      -       -       -       -
#populate per patient per column
for line in sys.stdin:
	line=line.strip().split()
	patient_id=line[0]
	cancer=line[1]
	#if new patient initialize
	if patient_id not in patient:
		patient[patient_id]=[]
		for x in range(cols):
			patient[patient_id].append(set())
	#for each thing update
	for x in range(cols):
		patient[patient_id][x].update(get_set(line[3+x]))
	#update genes_vs_tumors
	for x in range(cols):
		for gene in get_set(line[3+x]):
			if gene not in genes_vs_tumors[x]:
				genes_vs_tumors[x][gene]={}
			if cancer not in genes_vs_tumors[x][gene]:
				genes_vs_tumors[x][gene][cancer]=0
			genes_vs_tumors[x][gene][cancer]+=1
	#update tumors_vs_genes
	for x in range(cols):
		if cancer not in tumors_vs_genes[x]:
			tumors_vs_genes[x][cancer]={}
		for gene in get_set(line[3+x]):
			if gene not in tumors_vs_genes[x][cancer]:
				tumors_vs_genes[x][cancer][gene]=0
			tumors_vs_genes[x][cancer][gene]+=1
	# add to lookup
	if cancer not in cancers:
		cancers[cancer]=set()
	cancers[cancer].add(line[0])		
	if 'all' not in cancers:
		cancers['all']=set()
	cancers['all'].add(patient_id)		


#trim
trim=5
for x in range(cols):
	rt=set()
	for tumor in tumors_vs_genes[x]:
		v=set()
		for gene in tumors_vs_genes[x][tumor]:
			v.add(tumors_vs_genes[x][tumor][gene])
		if len(v)>2:
			v.remove(max(v))
			v.remove(max(v))
			etrim=max(trim,max(v))
		else:
			etrim=trim
		r=set()
		for gene in tumors_vs_genes[x][tumor]:
			if tumors_vs_genes[x][tumor][gene]<etrim:
				r.add(gene)
		for z in r:
			tumors_vs_genes[x][tumor].pop(z)
		if len(tumors_vs_genes[x][tumor])==0:
			rt.add(tumor)
	for t in rt:
		tumors_vs_genes[x].pop(t)

for x in range(cols):
	rg=set()
	for gene in genes_vs_tumors[x]:
		v=set()
		for tumor in genes_vs_tumors[x][gene]:
			v.add(genes_vs_tumors[x][gene][tumor])
		if len(v)>2:
			v.remove(max(v))
			v.remove(max(v))
			etrim=max(trim,max(v))
		else:
			etrim=trim
		r=set()
		for tumor in genes_vs_tumors[x][gene]:
			if genes_vs_tumors[x][gene][tumor]<etrim:
				r.add(tumor)
		for z in r:
			genes_vs_tumors[x][gene].pop(z)
		if len(genes_vs_tumors[x][gene])==0:
			rg.add(gene)
	for g in rg:
		genes_vs_tumors[x].pop(g)

#print tumors_vs_genes
for x in range(cols):
	print "###################"
	print "#TUMORS_VS_GENES_"+str(x)
	print "###################"
	#first get all genes 
	genes=set()
	for tumor in tumors_vs_genes[x]:
		genes.update(set(tumors_vs_genes[x][tumor].keys()))
	genes=list(genes)
	genes.sort()
	print "Tumor\t",
	print "\t".join(genes)
	for tumor in tumors_vs_genes[x]:
		print tumor+"\t",
		for gene in genes:
			if gene not in tumors_vs_genes[x][tumor]:
				print "0\t",
			else:
				print str(tumors_vs_genes[x][tumor][gene])+"\t",
		print ""

#print genes_vs_tumors
for x in range(cols):
	print "###################"
	print "#GENES_VS_TUMORS_"+str(x)
	print "###################"
	#first get all tumors
	tumors=set()
	for gene in genes_vs_tumors[x]:
		tumors.update(set(genes_vs_tumors[x][gene].keys()))
	tumors=list(tumors)
	tumors.sort()
	print "Gene\t",
	print "\t".join(tumors)
	for gene in genes_vs_tumors[x]:
		print gene+"\t",
		for tumor in tumors:
			if tumor not in genes_vs_tumors[x][gene]:
				print "0\t",
			else:
				print str(genes_vs_tumors[x][gene][tumor])+"\t",
		print ""

	
