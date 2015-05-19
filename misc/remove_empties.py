#!/usr/bin/python

import sys





def remove_empties(l,unscale=1):
	ll=[]
	u=set()
	for ia,ib,t,w in l:
		if w>0:
			w=int(w/unscale)
			ll.append((ia,ib,t,w))
			u.add(ia)
			u.add(ib)
	return ll,u
		

print >> sys.stderr, "#UNSCALE BY 1"

lines=sys.stdin.readlines()

lines=lines[-2:]


somatic=eval(lines[1])
somaticr,su=remove_empties(somatic,1)

genomic=eval(lines[0])
genomicr,gu=remove_empties(genomic,1)

for i in su:
	if i not in gu:
		print >> sys.stderr , "ERROR!"

print genomicr
print somaticr



