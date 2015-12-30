#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <zlib.h>
#include <omp.h>
using namespace std;

int lengths[26];

long reads_per_chr[26];
char ** fsta;

typedef struct __attribute__((__packed__)) struct_cov {
        unsigned short chr;
        unsigned int pos;
        unsigned short cov;
} struct_cov;

unsigned short get_chr(const char * s) {
	if (s[0]=='C' || s[0]=='c') {
		s+=3;
	}
	char tmp[10];
	int i=0; 
	for (i=0; i<10 && isalnum(s[i]); i++) {
		tmp[i]=s[i];
	}
	tmp[i]='\0';
	if (strlen(tmp)>0) {
		if (tmp[0]=='x' || tmp[0]=='X') {
			return 23;
		}
		if (tmp[0]=='y' || tmp[0]=='Y') {
			return 24;
		}
		if (tmp[0]=='m' || tmp[0]=='M') {
			return 25;
		}
		return atoi(tmp);
	} else {
		return 0;
	}
}

int gc(int chr, int pos, int size) {
	int gc=0;
	if (pos<=size || lengths[chr-1]<=(pos+size)) {
		return -1;
	}
	for (int i=0; i<size-1; i++) {
		char c = tolower(fsta[chr-1][pos+i-size/2]);
		switch (c) {
			case 'c':
			case 'g':
				gc++;
				break;
			case 'a':
			case 't':
				//meh
				break;
			case 'n':
				//meh
				return -1;
				break;
			default:
				cerr << "unexpected got char " << c << " | " << chr << ":" << pos << ":" << size << endl;
				return -1;
				//exit(1);
		}
	}
	return gc;
} 

long * read_cov(char * filename, bool normal, int bins) {
	long * gcbins = (long*)malloc(sizeof(long)*bins);
	if (gcbins==NULL) {
		cerr << "something bad" << endl;	
		exit(1);
	}

	for (int i=0; i<bins; i++) {
		gcbins[i]=0;
	}

	cerr << "Reading coverage from file " << filename << endl;

	gzFile fptr = gzopen(filename,"r");
	if (fptr==NULL) {
		fprintf(stderr, "Failed to open file %s\n",filename);
		exit(1);
	}
	
	//find the size of one entry
	size_t soe = sizeof(unsigned short)+sizeof(unsigned int)+sizeof(unsigned short);
	//size_t chunk = 30737418240L;
	size_t chunk = 1024*1024*1024L;
	size_t size_so_far = 0;
	char * buffer = (char*) malloc(chunk);
	if (buffer==NULL) {
		cerr << " FALLED TO MALLOC " << endl;
		exit(1);
	}


	for (int i=0; i<26; i++) {
		reads_per_chr[i]=0;
	}	
	
	//the real read loop
	while (!gzeof(fptr)) {
		size_t read = gzread(fptr,buffer+size_so_far,chunk);
		size_so_far+=read;
#ifdef DEBUG
		cerr << "Read so far " << size_so_far << endl;
#endif 
		if (read==chunk) {
			//cerr << "REALLOC" << endl;
			buffer=(char*)realloc(buffer,size_so_far+chunk);
			if (buffer==NULL) {
				cerr << " FALLED TO REALLOC " << endl;
				exit(1);
			}
		}
	}
	
	//lets get a buffer to fit the file	
	cout << "Done reading file  " << size_so_far <<  endl;

	size_t sz = size_so_far;	

	unsigned int entries = sz/soe;
	unsigned long total_coverage=0;

	cout << "Finding the average" << endl;


	
	int threads=32;
	omp_set_num_threads(threads); //omp_get_num_threads();

	long skips=0;

	long * t_reads_per_chr = (long*)malloc(threads*sizeof(long)*26);
	if (t_reads_per_chr==NULL) {
		cerr << "whoops" << endl;
		exit(1);
	}

	for (int i=0; i<threads*26; i++) {
		t_reads_per_chr[i]=0;
	}

	long * t_skips=(long*)malloc(threads*sizeof(long));
	if (t_skips==NULL) {
		cerr << "someting bad" << endl;
		exit(1);
	}
	for (int i=0; i<threads; i++) {
		t_skips[i]=0;
	}

	long * t_gcbins = (long*)malloc(threads*sizeof(long)*bins);
	if (t_gcbins==NULL) {
		cerr << "something bad" << endl;	
		exit(1);
	}
	for (int i=0; i<threads*bins; i++) {
		t_gcbins[i]=0;	
	}


	//setup threading counts for total coverage
	long * t_coverage = (long*)malloc(threads*sizeof(long)*bins);
	if (t_coverage==NULL) {
		cerr << "broke" << endl;
		exit(1);
	}
	for (int i=0; i<threads; i++) {
		t_coverage[i]=0;
	}

	cerr << "starting threading " << endl;

	#pragma omp parallel for schedule(static,1)
	for (unsigned int i=0; i<entries; i++) {
		int tidx=omp_get_thread_num();
		if (tidx>=threads) {
			cerr <<  "EMERGENCY STOP" << endl;
			exit(1);
		}
#ifdef DEBUG
		if (i%10000000==tidx) {
			cerr << tidx << " of " << omp_get_num_threads() << " : " << i << " / " << entries << endl;
		}
#endif
		if (omp_get_num_threads()!=threads) {
			cerr << "FAIL BAD" << endl;
			exit(1);
		}
		struct_cov * cov = (struct_cov*)(buffer+i*sizeof(struct_cov));
		if (cov->chr<26) {
			t_reads_per_chr[tidx*26+ cov->chr-1]+=cov->cov;
		}
		if (cov->chr>25) {
			cerr << "WHAT CHR" << cov->chr << endl;
		}
		int current_gc=gc(cov->chr,cov->pos,300);
		
		if (current_gc<0) {
			t_skips[tidx]++;
		} else {
			t_gcbins[tidx*bins + current_gc]+=cov->cov;
		}

		t_coverage[tidx]+=cov->cov;
		//total_coverage+=cov;
	}

	cerr << "merging " << endl;


	for (int i=0; i<threads; i++){ 
		total_coverage+=t_coverage[i];
	}	

	for (int i=0; i<threads; i++) {
		skips+=t_skips[i];
	}
	cerr << "SKIPPED: " << skips << endl;
	for (int i=0; i<threads; i++) {
		for (int j=0; j<bins; j++) {
			gcbins[j]+=t_gcbins[i*bins+j];
		}
		for (int j=0; j<26; j++) {
			reads_per_chr[j]+=t_reads_per_chr[i*26+j];
		}
	}


	cerr << "mereged" << endl;
	cout << "GC\t" ; 
	for (int i=0; i<bins; i++) {
		cout << gcbins[i] << "\t";
	}
	cout << endl;

	cout << "RPC\t" ;
	for (int i=0; i<26; i++) {
		cout << reads_per_chr[i] << "\t";
	}
	cout << endl;

	double average=((double)total_coverage)/entries;

	double sum=0;

	cout << "Average is " << average << " , now finding stddev" << endl;	
	for (unsigned int i=0; i<entries; i++) {
		char* base = buffer+i*soe;
		unsigned short chr=*((unsigned short *)base);
		if (chr==0) {
			continue;
		}
		base+=sizeof(unsigned short);
		unsigned int coord=*((unsigned int *)base);
		base+=sizeof(unsigned int);
		unsigned short cov=*((unsigned short *)base);


		sum+=(cov-average)*(cov-average)/entries;
	}


	double stddev=sqrt(sum);
	
	cout << "Standard deviation is " << stddev << endl;

	


	cerr << "total: " << total_coverage << endl;
	free(buffer);

	return gcbins;
}


char ** read_fasta(char * filename) {
	fsta = (char**)malloc(sizeof(char*)*26);
	if (fsta==NULL) {
		cerr << "NOT GOOD" << endl;
		exit(1);
	}

	gzFile  f = gzopen(filename,"r");
	if (f==NULL) {
		cerr << " NOT GOOD 2 " << endl;
	}


	char * ref = (char*)malloc(sizeof(char)*5000000000L);
	if (ref==NULL) {
		fprintf(stderr, "ERROR\n");
		exit(2);
	}

	cout << "-Reading reference" << endl;
	size_t ret=0;
	int inc=0;
	inc = gzread(f,ref,20000000);
	while (inc>0) {
		ret+=inc;
		inc = gzread(f,ref+ret,20000000);
	}
	//size_t ret = gzread(f,ref,20000000000L);
	cout << "Read " << ret << " from ref" << endl;
	cout << "+Done reading reference" << endl;


	char * buffer = (char*)malloc(sizeof(char)*1000000000L);
	size_t i =0; 
	size_t j =0;
	buffer[j]='\0';
	int ichr=0;
	bool in_header=false;
	for (i=0; i<ret ; i++) {
		if (in_header || ref[i]=='\n') {
			if (in_header && ref[i]=='\n') {
				in_header=false;
			}
			continue;
		} else if (ref[i]=='>') {
			cerr << "-Found header line for chr " << get_chr(ref+i+1) << endl;
			//copy out the old chromosome
			if (ichr!=0) {
				fsta[ichr-1]=(char*)malloc(sizeof(char)*(j+1));
				if (fsta[ichr-1]==NULL) {
					cerr << "MAJOR ERROR" << endl;
					exit(1);
				}
				memcpy(fsta[ichr-1],buffer,j*sizeof(char));
				fsta[ichr-1][j]='\0';
				cout << "+Processed chr " << ichr << " at index " << i << endl;
				char  tmp[100];
				memcpy(tmp,fsta[ichr-1]+1000000,99);
				tmp[99]='\0';
				cout << tmp << "...";
				memcpy(tmp,&fsta[ichr-1][j-100],99);
				tmp[99]='\0';
				cout << tmp << endl;
				lengths[ichr-1]=strlen(fsta[ichr-1]);
				for (int x=0; x< lengths[ichr-1]; x++) {
					switch(tolower(fsta[ichr-1][x])) {
						case 'a':
						case 'c':
						case 't':
						case 'g':
						case 'n':
							break;
						default:
							cerr << "unknown char |" << fsta[ichr-1][x] << "| at " << x << endl;
					}
				}		
				j=0;
				
			}
			//start working on next chromsome
			ichr=get_chr(ref+i+1);
			in_header=true;
		} else {
			buffer[j++]=ref[i];
		}
	}
	if (ichr!=0) {
		fsta[ichr-1]=(char*)malloc(sizeof(char)*(j+1));
		if (fsta[ichr-1]==NULL) {
			cerr << "MAJOR ERROR" << endl;
			exit(1);
		}
		memcpy(fsta[ichr-1],buffer,j*sizeof(char));
		fsta[ichr-1][j]='\0';
		lengths[ichr-1]=strlen(fsta[ichr-1]);
				for (int x=0; x< lengths[ichr-1]; x++) {
					switch(tolower(fsta[ichr-1][x])) {
						case 'a':
						case 'c':
						case 't':
						case 'g':
						case 'n':
							break;
						default:
							cerr << "unknown char |" << fsta[ichr-1][x] << "| at " << x << endl;
					}
				}		
		cout << "+Processed chr " << ichr << " at index " << i << endl;
	}
	cout << "Finished reference processing" << endl;
	for (int i=0; i<24; i++) {
		cerr << "chr "<< i+1 << lengths[i] << endl;
	}
	return fsta;	
}





int main (int argc, char ** argv) {
	if (argc!=3) {
		cerr<<argv[0]<<" in_file ref"<<endl;
		return 0;
	}

	char * coverage_filename=argv[1];
	char * fasta_filename=argv[2];

	read_fasta(fasta_filename);
	read_cov(coverage_filename,false,300);

}

//chr1    10000   N       3       AA^:A   D@+
