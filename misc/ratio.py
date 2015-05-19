#!/usr/bin/python

import sys

if len(sys.argv)!=1:
	print "crash"
	sys.exit(1)


def count(s,z,t):
	x=s.replace(z.upper(),'').replace(z.lower(),'')
	if t:
		x=x.replace('.','').replace(',','')
	return len(s)-len(x)

for line in sys.stdin:
	#print line
	line = line.split()
	ref = line[2].lower()
	a = count(line[4],'a',ref=='a')
	c = count(line[4],'c',ref=='c')
	g = count(line[4],'g',ref=='g')
	t = count(line[4],'t',ref=='t')
	rc=0
	if ref=='a':
		rc=a
	elif ref=='c':
		rc=c
	elif ref=='g':
		rc=g
	elif ref=='t':
		rc=t	
	ratio = float(rc)/(a+c+g+t+1)
	print "\t".join([line[0],line[1],"%0.3f" % ratio,str(a),str(c),str(g),str(t)])
	#chr1    65872   T       56      ,gg..GG.GGg,G,G,,,g,.Ggg,g,G,,,Gggg,,gg,G,G,GGG,GGGG,g,^!g
