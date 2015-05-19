import sys


next_binary=False
next_integer=False
for line in sys.stdin:
	if line=="Minimize\n":
		continue
	elif line[:4]==" obj":
		line=line.replace(" obj","min")
		print line.strip()+";"
	elif line=="Subject To\n":
		continue
	elif line=="Bounds\n":
		continue
	elif line=="End\n":
		continue
	elif line=="Binary\n":
		next_binary=True
	elif next_binary:
		next_binary=False
		print "bin "+", ".join(line.strip().split())+";" 
	elif line=="General\n":
		next_integer=True
	elif next_integer:
		next_integer=False
		#print "int "+", ".join(line.strip().split())+";" 
	elif line[0]=='\\':
		print "/*"+line.strip()+"*/"
	else:
		print line.strip()+";"
