import sys
import xml.etree.ElementTree as ET

if len(sys.argv)!=2:
	print sys.argv[0], "xml_filename"
	sys.exit(1)


xml_filename=sys.argv[1]

tree = ET.parse(xml_filename)
root = tree.getroot()


dna_samples={}



for c in root.findall("./Result"):
	r=c.findall("./analyte_code")
	r2=c.findall("./center_name")
	r3=c.findall("./sample_type")
	if len(r)==1 and len(r2)==1 and len(r3)==1:
		sample_type=int(r3[0].text)
		center_name=r2[0].text
		if center_name not in dna_samples:
			dna_samples[center_name]={}
		if r[0].text=="D":
			dna_samples[center_name][sample_type]=c
	else:
		print r
		continue

if len(dna_samples)==0:
	sys.exit(1)

largest=0
largest_pair=[]


#trim by type
for center in dna_samples:
	samples=dna_samples[center]
	if len(samples)>2:
		#get the smallest id lower then 10
		#get the smallest id larger or equal to 10
		tumor_id=100
		normal_id=100
		for k in samples.keys():
			if k<10 and k<tumor_id:
				tumor_id=k
			if k>=10 and k<normal_id:
				normal_id=k
		samples={tumor_id:samples[tumor_id],normal_id:samples[normal_id]}
	size=0
	for sid in samples:
		y=samples[sid]
		size+=sum(map(lambda x : int(x.text), y.findall('./files/file/filesize')))
	if size>largest:
		largest=size
		largest_pair=samples

if len(largest_pair)==2:
	for sid in largest_pair:
		print largest_pair[sid].findall('./analysis_id')[0].text

	
