#/usr/bin/python


import sys


if len(sys.argv)!=3:
	print "%s p f - result is p-f" % sys.argv[0]
	sys.exit(1)


pfn=sys.argv[1]
flfn=sys.argv[2]


def read_flow(fn):
	d={}
	f=open(fn,'r')
	for line in f:
		if line[0]=='f':
			line=line.split()
			s=int(line[1])
			t=int(line[2])
			l=int(line[3])
			if not s in d:
				d[s]={}
			if t not in d[s]:
				d[s][t]=0
			d[s][t]+=l		
	f.close()
	return d

def read_problem_file(fn,fl):
	outb=[]
	f=open(fn,'r')
	b=[]
	arcs=0
	for line in f:
		if line[0]=='n':
			if len(b)>0:
				outb.append("\n".join(b))
				del b[:]
			outb.append(line.strip())
		elif line[0]=='c':
			if len(b)>0 and b[-1][0]=='c':
				del b[:]
			b.append(line.strip())
		elif line[0]=='a':
			line_=line.split()
			s=int(line_[1])
			t=int(line_[2])
			l=int(line_[4])
			if s in fl and t in fl[s] and fl[s][t]>0:
				#need to skip this entry
				k=min(fl[s][t],l)
				l-=k
				fl[s][t]-=k
			if l>0:
				if len(b)>0:
					outb.append("\n".join(b))
				line_[4]=str(l)
				outb.append("\t".join(line_).strip())
				arcs+=1
			if len(b)>0:
				del b[:]
		elif line[0]=='p':
			if len(b)>0:
				outb.append("\n".join(b))
				del b[:]
			outb.append(line.strip())
	f.close()
	for line in outb:
		if line[0]=='p':
			line=line.split()
			line[3]=str(arcs)
			print "\t".join(line)
		else:
			print line

fl=read_flow(flfn)
ks=fl.keys()
for s in ks:
	kt=fl[s].keys()
	for t in kt:
		if not t in fl:
			fl[t]={}
		if not s in fl[t]:
			fl[t][s]=0
		fl[t][s]=max(fl[t][s],fl[s][t])
		fl[s][t]=fl[t][s]

read_problem_file(pfn,fl)

