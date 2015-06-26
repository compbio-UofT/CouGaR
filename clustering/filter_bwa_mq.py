#!/usr/bin/python


import sys

if len(sys.argv)!=2:
	print "%s min_mapq" % sys.argv[0]

min_mq=int(sys.argv[1])

for line in sys.stdin:
	if line[0]=='@':
		print line,
	elif line.find("MQ:i:")>0:
		mq=int(line.split("MQ:i:")[1].split()[0])
		my_mq=int(line.split('\t')[4])
		if mq>=min_mq and my_mq>=min_mq:
			print line,
	else:
		print line,

		
		
