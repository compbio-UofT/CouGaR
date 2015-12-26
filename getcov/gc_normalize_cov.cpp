#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <zlib.h>
#include <omp.h>
#include <assert.h>
using namespace std;


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MAX_MEDIAN 1000

int lengths[26];

long reads_per_chr_tumor[26];
long reads_per_chr_normal[26];
int normal_median; 
int tumor_median;
long total_tumor; 
long total_normal;
long * normal_gcbins;
long * tumor_gcbins;
char ** fsta;
int ** gcs;

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


void correct_tumor_coverage(char * buffer, unsigned int entries) {
	#pragma omp parallel for 
	for (unsigned int i=0; i<entries; i++) {
		int tidx=omp_get_thread_num();
		if (i%100000000==tidx) {
			cerr << tidx << " of " << omp_get_num_threads() << " : " << i << " / " << entries << endl;
		}
		struct_cov * cov = (struct_cov*) (buffer+i*sizeof(struct_cov));
		if (cov->chr>25) {
			cerr << "WHAT CHR" << cov->chr << endl;
		}
		//int current_gc=gc(cov->chr,cov->pos,300);
		int current_gc = gcs[cov->chr-1][cov->pos];
		assert( current_gc>=-1);	
		if (current_gc>=0) {
			double gc_normal = ((double)normal_gcbins[current_gc])/total_normal;
			double gc_tumor = ((double)tumor_gcbins[current_gc])/total_tumor;
			cov->cov=(cov->cov)*(gc_tumor/gc_normal)*(tumor_median/normal_median); // correct for GC BIN and coverage
		}
	}
}

long * read_cov(char * filename, int bins, long * reads_per_chr, char ** ret_buffer, unsigned int * ret_entries, int * median) {
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
		cerr << "Read so far " << size_so_far << endl;
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
	assert(soe==sizeof(struct_cov));
	unsigned long total_coverage=0;

	cout << "Finding the average from " << entries << endl;


	
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


	int * t_median_counts = (int*)malloc(threads*sizeof(int)*MAX_MEDIAN);
	if (t_median_counts==NULL) {
		cerr << "FAILED TO GET COUNTS FOR MEDIAN" << endl;
		exit(1);
	}
	memset(t_median_counts,0,sizeof(int)*MAX_MEDIAN*threads);

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

	#pragma omp parallel for 
	for (unsigned int i=0; i<entries; i++) {
		int tidx=omp_get_thread_num();
		/*if (tidx>=threads) {
			cerr <<  "EMERGENCY STOP" << endl;
			exit(1);
		}*/
		if (i%100000000==tidx) {
			cerr << tidx << " of " << omp_get_num_threads() << " : " << i << " / " << entries << endl;
		}
		//if (omp_get_num_threads()!=threads) {
		//	cerr << "FAIL BAD" << endl;
		//	exit(1);
		//}
		/*char* base = buffer+i*soe;
		unsigned short chr=(*((unsigned short *)base));
		if (chr==0) {
			t_skips++;
			continue;
		}
		base+=sizeof(unsigned short);
		unsigned int coord=*((unsigned int *)base);
		base+=sizeof(unsigned int);
		unsigned short cov=*((unsigned short *)base);*/
		struct_cov * cov = (struct_cov*) (buffer+i*sizeof(struct_cov));
		if (cov->chr<26) {
			t_reads_per_chr[tidx*26+ cov->chr-1]+=cov->cov;
		}
		if (cov->chr>25) {
			cerr << "WHAT CHR" << cov->chr << endl;
		}

		int median_cov = MIN(MAX_MEDIAN-1,cov->cov);
		t_median_counts[tidx*MAX_MEDIAN+median_cov]++;
		//int current_gc=gc(cov->chr,cov->pos,300);
		int current_gc = gcs[cov->chr-1][cov->pos];
		assert(current_gc<=bins && current_gc>=-1);	
		
		if (current_gc<0) {
			t_skips[tidx]++;
		} else {
			t_gcbins[tidx*bins + current_gc]+=cov->cov;
		}

		t_coverage[tidx]+=cov->cov;
		//total_coverage+=cov;
	}

	cerr << "merging " << endl;

	for (int i=1; i<threads; i++) { 
		for (int j=0; j<MAX_MEDIAN; j++) {
			t_median_counts[j]+=t_median_counts[i*MAX_MEDIAN+j];
		}
	}	
	int max_median_idx = 0;
	for (int j=0; j<MAX_MEDIAN; j++) {
		if ( t_median_counts[max_median_idx] < t_median_counts[j] ) {
			max_median_idx = j;
		}
	}
	cerr << "MEDIAN " << max_median_idx << endl;

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
	//free(buffer)
	*ret_buffer =buffer;
	*ret_entries = entries;

	free(t_median_counts);
	return gcbins;
}


char ** read_fasta(char * filename, int gc_size) {
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
		cout << "chr"<< i+1  << " " << lengths[i] << endl;
	}

	//compute the gcs
	gcs=(int**)malloc(sizeof(int*)*25);
	if (gcs==NULL) {
			cerr << "failed malloc for gcs" << endl;
			exit(1);
	}

	//#pragma omp parallel for 
	for (int i=0; i<25; i++) {
		gcs[i] = (int*)(malloc(sizeof(int)*lengths[i]));
		if (gcs[i]==NULL) {
			cerr << "failed malloc for gcs" << endl;
			exit(1);
		}
		
		int gc=0;
		int at=0;
		int ns=0;
		for (int j=0; j<lengths[i]; j++) {
			int idx = j-gc_size/2;
			char c = tolower(fsta[i][j]);
			switch (c) {
				case 'c':
				case 'g':
					gc++;
					break;
				case 'a':
				case 't':
					at++;
					break;
				case 'n':
					ns++;
					break;
				default:
					cerr << "unexpected got char " << c << endl;
					exit(1);
			}
			int ridx = j-gc_size;
			if (ridx>=0) {
				char c = tolower(fsta[i][ridx]);
				switch (c) {
					case 'c':
					case 'g':
						gc--;
						break;
					case 'a':
					case 't':
						at--;
						break;
					case 'n':
						ns--;
						break;
					default:
						cerr << "unexpected got char " << c << endl;
						exit(1);
				}
			}
			if (j<gc_size) {
				if (idx>=0) {
					gcs[i][idx]=-1;
				}
			} else {
				if (ns>0) {
					gcs[i][idx]=-1;
				} else {
					gcs[i][idx]=gc;
				}
			}
		}
		for (int idx=lengths[i]-gc_size/2; idx<lengths[i]; idx++) {
			gcs[i][idx]=-1;
		}
	}
	return fsta;	
}





int main (int argc, char ** argv) {
	if (argc!=5) {
		cerr<<argv[0]<<" in_file_normal in_file_tumor ref output_filename"<<endl;
		return 0;
	}

	char * normal_coverage_filename=argv[1];
	char * tumor_coverage_filename=argv[2];
	char * fasta_filename=argv[3];
	char * output_filename = argv[4];

	int bins = 300;
	read_fasta(fasta_filename,bins);
	char * buffer = NULL;
 	unsigned int entries;
	normal_gcbins = read_cov(normal_coverage_filename,bins,reads_per_chr_normal, &buffer, &entries, &normal_median);
	free(buffer);
	tumor_gcbins = read_cov(tumor_coverage_filename,bins,reads_per_chr_tumor, &buffer, &entries, &tumor_median);

	total_normal = 0;
	total_tumor = 0;
	for (int i=0; i<bins; i++) {
		total_normal+=normal_gcbins[i];
		total_tumor+=tumor_gcbins[i];
	}

	//lets normalize the tumor sample to be like the normal?	
	for (int i=0; i<bins; i++) {
		double gc_normal = ((double)normal_gcbins[i])/total_normal;
		double gc_tumor = ((double)tumor_gcbins[i])/total_tumor;
		cout << i << " " <<  gc_normal << " " << gc_tumor << endl;
	}

	
	cerr << "writting file" << endl;
	gzFile *fi = (gzFile *)gzopen(output_filename,"wb");
        if (fi==NULL) {
                cerr<<"something terrible";
                return 1;
        }
	unsigned int i=0;
	int write_at_once = 1024*1024*32;
	while (i<entries) {
		unsigned int left = entries - i;
		int to_write=write_at_once;
		if (to_write>left) {
			to_write=left;
		}
		gzwrite(fi,buffer+sizeof(struct_cov)*i,sizeof(struct_cov)*to_write);
		cerr << i << " " << to_write << " " << ((double)i/entries) << endl;
		i+=to_write;
	}

	gzclose(fi);	
	
	
}

//chr1    10000   N       3       AA^:A   D@+
