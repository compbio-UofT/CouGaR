import sys
import gzip


if len(sys.argv)!=3:
	print "%s sandborn walk" % sys.argv[0]
	sys.exit(1)


def overlap(sbr,wr):
	if sbr[0]!=wr[0]:
		return False
	#check when sbr is subsumed by wr
	if sbr[1]>(wr[1]-100) and sbr[2]<(wr[2]+100):
		return True

def bp_overlap(sbr,wr):
	mxs=max(sbr[1],wr[1])
	mne=min(sbr[2],wr[2])
	if mxs<mne:
		return mne-mxs
	return 0 

def intersect(sbr,wr):
	mxs=max(sbr[1],wr[1])
	mne=min(sbr[2],wr[2])
	if mxs<mne:
		return True
	return False

def to_coord(x):
	s=str(x)
	ns=""
	for z in range(len(s)):
		if z>0 and z%3==2:
			ns=ns+","
		ns=ns+s[z]
	return ns
		
	

def to_bp(x):
	if x>=100000:
		return "%0.2fMb" % (x/float(1000000))
	if x>=1000:
		return "%0.2fKb" % (x/float(1000))
	return "%db" % x
	

def resolve(sb,w):
	for wr in w:
		o=[]
		o.append("%s:%s+%s\t" % (wr[0],to_coord(wr[1]),to_bp(wr[2]-wr[1])))
		intersects=[]
		for sbr in sb:
			if intersect(sbr,wr):
				intersects.append(sbr)
		#lets walk
		c=[wr[0],wr[1],wr[2]]
		for sbr in intersects:
			if sbr[1]>c[1]:
				#print str(sbr[1]-c[1]) +"bp,",
				o.append(to_bp(sbr[1]-c[1]))
				c[1]=sbr[1]
			if sbr[1]<=c[1]:
				fc=(100*bp_overlap(sbr,wr)/float(sbr[2]-sbr[1]))
				o.append(sbr[3]+("(%0.1f)" % fc))
				c[1]=sbr[2]
		if c[1]<c[2]:
			#print str(c[2]-c[1]) +"bp,",
			o.append(to_bp(c[2]-c[1]))
		print o[0],
		if len(o)>1:
			print ",".join(o[1:])

				
				
		


sandborn_filename=sys.argv[1]
walk_filename=sys.argv[2]

sandborn_regions=[]

sandborn_file=open(sandborn_filename)
for line in sandborn_file:
	line=line.strip().split()
	chr=line[0]
	s=int(line[1])
	e=int(line[2])
	l=line[3]
	if s>e:
		print "ERROR"
		sys.exit(1)
	sandborn_regions.append((chr,s,e,l))

sandborn_regions.sort()

walk_regions=[]

walk_file=open(walk_filename)
#(12, 58168188) (12, 58173832) *
for line in walk_file:
	line=line.replace(")",'').replace("(",'').replace(",",'').replace(" *",'').strip().split()
	chr="chr"+line[0]
	s=int(line[1])
	e=int(line[3])
	if e<s:
		t=e
		e=s
		s=t
	if s>e:
		print "ERROR"
		sys.exit(1)
	walk_regions.append((chr,s,e))




resolve(sandborn_regions,walk_regions)
