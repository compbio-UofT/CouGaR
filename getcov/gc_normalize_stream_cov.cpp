#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <zlib.h>
#include <omp.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
using namespace std;


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MAX_MEDIAN 1000

int threads = 32;
int bins = 300;
int lengths[26];
unsigned long genome_length =0; 
int normal_median; 
int tumor_median;
char ** fsta;
int ** gcs;

typedef struct struct_stats {
	long * gcbins;
	long * reads_per_chr;
	double mean;
	double stddev;
	int median;	
	int skips;
	unsigned long total_reads;
	double mean_bar;
	double stddev_bar;
} struct_stats;

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


void correct_coverage(char * buffer, unsigned int entries, bool tumor_b, struct_stats * normal_stats, struct_stats *tumor_stats) {
	cerr << "Starting correction for " << entries << " entries" << endl;
	//double median_ratio = ((double)normal_median/(double)tumor_median);
	double median_ratio = 1.0;
	if (tumor_b) {
		median_ratio = normal_stats->mean/tumor_stats->mean;	
	} else {
		median_ratio = tumor_stats->mean/normal_stats->mean;	
	}
	if (median_ratio>1.0) {
		median_ratio=1.0;
	}
	srand(time(NULL));
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
			double gc_ratio = 1.0;
			if (tumor_b) {
				double gc_normal = ((double)normal_stats->gcbins[current_gc])/normal_stats->total_reads;
				double gc_tumor = ((double)tumor_stats->gcbins[current_gc])/tumor_stats->total_reads;
				gc_ratio = gc_normal/((1-gc_normal)*gc_tumor);
			} else {
				double gc_normal = ((double)normal_stats->gcbins[current_gc])/normal_stats->total_reads;
				double gc_tumor = ((double)tumor_stats->gcbins[current_gc])/tumor_stats->total_reads;
				gc_ratio = gc_tumor/((1-gc_tumor)*gc_normal);
			}
			if (gc_ratio>1) {
				gc_ratio=1.0;
			}
			//cerr << cov->cov << " " << gc_ratio << " " << median_ratio << " " << cov->cov*gc_ratio*median_ratio << endl;
			double new_cov = cov->cov*gc_ratio*median_ratio;
			double X=((double)rand()/(double)RAND_MAX);
			if (X<(new_cov - (int)new_cov) ) {
				new_cov++;
			}
			cov->cov=(int)new_cov; // correct for GC BIN and coverage

		}
	}
	cerr << "Correction complete" << endl;
}


struct_stats * new_stats() {
	struct_stats * ss = (struct_stats*)malloc(sizeof(struct_stats));
	if (ss==NULL) {
		cerr << "Failed to allocate new stats" << endl;
		exit(1);
	}
	memset(ss,0,sizeof(struct_stats));

	ss->gcbins = (long*)malloc(sizeof(long)*bins);
	if (ss->gcbins==NULL) {
		cerr << "failed to allocate new stats" << endl;
		exit(1);
	}
	memset(ss->gcbins,0,sizeof(long)*bins);
	
	ss->reads_per_chr = (long*)malloc(sizeof(long)*26);
	if (ss->reads_per_chr == NULL) {
		cerr << "failed to allocate new stast x2" << endl;
		exit(1);
	}
	memset(ss->reads_per_chr,0,sizeof(long)*26);

	//probably dont need this after top memset, but why not
	ss->mean=0.0;
	ss->stddev=0.0;
	ss->median=0;	
	ss->skips=0;
	ss->total_reads=0;
	ss->mean_bar=0;
	ss->stddev_bar=0;
	return ss;
}

void free_stats(struct_stats *ss) {
	free(ss->gcbins);
	free(ss->reads_per_chr);
	free(ss);
}

struct_stats *  find_stats(char * buffer, unsigned int entries) {
	//cerr << "Starting to find stats" << endl;
	struct_stats * ss = new_stats();

	//thread local storage
	long * t_reads_per_chr = (long*)malloc(threads*sizeof(long)*26);
	if (t_reads_per_chr==NULL) {
		cerr << "whoops" << endl;
		exit(1);
	}
	memset(t_reads_per_chr, 0, sizeof(long)*threads*26);

	long * t_skips=(long*)malloc(threads*sizeof(long));
	if (t_skips==NULL) {
		cerr << "someting bad" << endl;
		exit(1);
	}
	memset(t_skips,0,sizeof(long)*threads);

	long * t_gcbins = (long*)malloc(threads*sizeof(long)*bins);
	if (t_gcbins==NULL) {
		cerr << "something bad" << endl;	
		exit(1);
	}
	memset(t_gcbins,0,sizeof(long)*threads*bins);

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
	memset(t_coverage,0,sizeof(long)*threads*bins);

	#pragma omp parallel for 
	for (unsigned int i=0; i<entries; i++) {
		int tidx=omp_get_thread_num();
		struct_cov * cov = (struct_cov*) (buffer+i*sizeof(struct_cov));
		if (cov->chr<26) {
			t_reads_per_chr[tidx*26+ cov->chr-1]+=cov->cov;
		}
		if (cov->chr>25) {
			cerr << "WHAT CHR" << cov->chr << endl;
		}
	
		//find median stats
		int median_cov = MIN(MAX_MEDIAN-1,cov->cov);
		t_median_counts[tidx*MAX_MEDIAN+median_cov]++;

		//find the current gcbin
		int current_gc = gcs[cov->chr-1][cov->pos];
		assert(current_gc<=bins && current_gc>=-1);	
		
		if (current_gc<0) {
			t_skips[tidx]++;
		} else {
			t_gcbins[tidx*bins + current_gc]+=cov->cov;
		}

		t_coverage[tidx]+=cov->cov;
	}

	//merge the results
	for (int i=0; i<threads; i++){ 
		ss->total_reads+=t_coverage[i];
		ss->skips+=t_skips[i];
		if (i>0) {
			for (int j=0; j<MAX_MEDIAN; j++) {
				t_median_counts[j]+=t_median_counts[i*MAX_MEDIAN+j];
			}
		} else {
			cerr << genome_length << " " << entries << endl;
			if (entries<genome_length) {
				t_median_counts[0]=genome_length - entries;
			}
		}
		for (int j=0; j<bins; j++) {
			ss->gcbins[j]+=t_gcbins[i*bins+j];
		}
		for (int j=0; j<26; j++) {
			ss->reads_per_chr[j]+=t_reads_per_chr[i*26+j];
		}
	}

	//mean
	//ss->mean=((double)ss->total_reads)/entries;
	ss->mean=((double)ss->total_reads)/genome_length; //correct for full genome length

	//median
	int max_median_idx = 0;
	unsigned long reads_so_far=0;
	for (int j=0; j<MAX_MEDIAN; j++) {
		if (j<30) {
			cerr << "BIN DIST X " << j << " " << t_median_counts[j] << endl;
		}
		reads_so_far+=t_median_counts[j];
		if ( (entries/2)<reads_so_far ) {
			ss->median=j;
			break;
		}
	}

	//stddev
	double sum=0;
	//cout << "Average is " << average << " , now finding stddev" << endl;	
	for (unsigned int i=0; i<entries; i++) {
		struct_cov * cov = (struct_cov*) (buffer+i*sizeof(struct_cov));
		sum+=(cov->cov-ss->mean)*(cov->cov-ss->mean)/entries;
	}
	ss->stddev=sqrt(sum);


	//
	//
	//
	//Second pass
	//
	//

	//do a second pass for mean_bar and stddev_bar
	memset(t_coverage,0,sizeof(long)*threads*bins);
	#pragma omp parallel for 
	for (unsigned int i=0; i<entries; i++) {
		int tidx=omp_get_thread_num();
		struct_cov * cov = (struct_cov*) (buffer+i*sizeof(struct_cov));
		if (cov->chr>25) {
			cerr << "WHAT CHR" << cov->chr << endl;
		}
		if ( (cov->cov-ss->mean)/ss->stddev<3 )  {
			t_coverage[tidx]+=cov->cov;
		}
	}

	//merge the results
	unsigned long total_reads = 0;
	for (int i=0; i<threads; i++) { 
		total_reads+=t_coverage[i];
	}

	//mean_bar
	ss->mean_bar=((double)total_reads)/entries;

	//stddev_bar
	sum=0;
	for (unsigned int i=0; i<entries; i++) {
		struct_cov * cov = (struct_cov*) (buffer+i*sizeof(struct_cov));
		sum+=(cov->cov-ss->mean_bar)*(cov->cov-ss->mean_bar)/entries;
	}
	ss->stddev_bar=sqrt(sum);

	//free temporary thread stuff
	free(t_reads_per_chr);
	free(t_skips);
	free(t_gcbins);
	free(t_median_counts);
	free(t_coverage);
	//cerr << "Done find stats" << endl;
	return ss;
}


void print_stats(struct_stats * ss ) {
	//print out the GC
	cout << "GC\t" ; 
	for (int i=0; i<bins; i++) {
		cout << ss->gcbins[i] << "\t";
	}
	cout << endl;

	//print out the RPC
	cout << "RPC\t" ;
	for (int i=0; i<26; i++) {
		cout << ss->reads_per_chr[i] << "\t";
	}
	cout << endl;

	cout << "Mean\t" <<  ss->mean << "\tStd\t" << ss->stddev << endl;
	cout << "Mean Bar\t" <<  ss->mean_bar << "\tStd_bar\t" << ss->stddev_bar << endl;
	cout << "Median\t" << ss->median << endl;
}


struct_stats * stream_read_cov(char * filename, int bins, bool correct, struct_stats * normal_stats, struct_stats * tumor_stats, char * output_filename, bool tumor_b) {
	//cerr << "Starting to find stats" << endl;
	struct_stats * ss = new_stats();

	//thread local storage
	long * t_reads_per_chr = (long*)malloc(threads*sizeof(long)*26);
	if (t_reads_per_chr==NULL) {
		cerr << "whoops" << endl;
		exit(1);
	}
	memset(t_reads_per_chr, 0, sizeof(long)*threads*26);

	long * t_skips=(long*)malloc(threads*sizeof(long));
	if (t_skips==NULL) {
		cerr << "someting bad" << endl;
		exit(1);
	}
	memset(t_skips,0,sizeof(long)*threads);

	long * t_gcbins = (long*)malloc(threads*sizeof(long)*bins);
	if (t_gcbins==NULL) {
		cerr << "something bad" << endl;	
		exit(1);
	}
	memset(t_gcbins,0,sizeof(long)*threads*bins);

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
	memset(t_coverage,0,sizeof(long)*threads*bins);



	//cerr << "Reading coverage from file " << filename << endl;

	gzFile fptr = gzopen(filename,"r");
	if (fptr==NULL) {
		fprintf(stderr, "Failed to open file %s\n",filename);
		exit(1);
	}
	
	//find the size of one entry
	//size_t chunk = 30737418240L;
	size_t chunk = 512*1024*1024L;
	char * buffer = (char*) malloc(chunk);
	if (buffer==NULL) {
		cerr << " FALLED TO MALLOC " << endl;
		exit(1);
	}

	double median_ratio = 1.0;
	if (correct) {
		if (tumor_b) {
			median_ratio = normal_stats->mean/tumor_stats->mean;	
		} else {
			median_ratio = tumor_stats->mean/normal_stats->mean;	
		}
		if (median_ratio>1.0) {
			median_ratio=1.0;
		}
	}
	srand(time(NULL));



	gzFile *fi = NULL;
	if (output_filename!=NULL) {
		fi = (gzFile *)gzopen(output_filename,"wb");
		if (fi==NULL) {
			cerr<<"something terrible";
			exit(1);
		}
	}



	//the real read loop
	unsigned long total_entries  = 0;
	while (!gzeof(fptr)) {
		size_t read = gzread(fptr,buffer,chunk);
		if (read%sizeof(struct_cov)!=0) {
			cerr << "FAILED TO READ FILE PROPERLY..." << endl;
			exit(1);
		}
		unsigned long entries = read/sizeof(struct_cov);
		total_entries+=entries;
		#pragma omp parallel for 
		for (unsigned int i=0; i<entries; i++) {
			int tidx=omp_get_thread_num();
			struct_cov * cov = (struct_cov*) (buffer+i*sizeof(struct_cov));

			if (cov->chr<26) {
				t_reads_per_chr[tidx*26+ cov->chr-1]+=cov->cov;
			}
			if (cov->chr>25) {
				cerr << "WHAT CHR" << cov->chr << endl;
				continue;
			}
			if (cov->pos>=lengths[cov->chr-1]) {
				cerr << "falling off the chromosome... " << cov->chr << " " << cov->pos << " END is " << lengths[cov->chr-1] << endl;
				continue;
			}
	
			if (correct) {	
				//int current_gc=gc(cov->chr,cov->pos,300);
				int current_gc = gcs[cov->chr-1][cov->pos];
				assert( current_gc>=-1);	
				if (current_gc>=0) {
					double gc_ratio = 1.0;
					if (tumor_b) {
						double gc_normal = ((double)normal_stats->gcbins[current_gc])/normal_stats->total_reads;
						double gc_tumor = ((double)tumor_stats->gcbins[current_gc])/tumor_stats->total_reads;
						gc_ratio = gc_normal/((1-gc_normal)*gc_tumor);
					} else {
						double gc_normal = ((double)normal_stats->gcbins[current_gc])/normal_stats->total_reads;
						double gc_tumor = ((double)tumor_stats->gcbins[current_gc])/tumor_stats->total_reads;
						gc_ratio = gc_tumor/((1-gc_tumor)*gc_normal);
					}
					if (gc_ratio>1) {
						gc_ratio=1.0;
					}
					//cerr << cov->cov << " " << gc_ratio << " " << median_ratio << " " << cov->cov*gc_ratio*median_ratio << endl;
					double new_cov = cov->cov*gc_ratio*median_ratio;
					//double X=((double)rand()/(double)RAND_MAX);
					//if (X<(new_cov - (int)new_cov) ) {
					//	new_cov++;
					//}
					cov->cov=(int)new_cov; // correct for GC BIN and coverage

				}
			}
		
			//find median stats
			int median_cov = MIN(MAX_MEDIAN-1,cov->cov);
			t_median_counts[tidx*MAX_MEDIAN+median_cov]++;

			//find the current gcbin
			int current_gc = gcs[cov->chr-1][cov->pos];
			assert(current_gc<=bins && current_gc>=-1);	
			
			if (current_gc<0) {
				t_skips[tidx]++;
			} else {
				t_gcbins[tidx*bins + current_gc]+=cov->cov;
			}

			t_coverage[tidx]+=cov->cov;
		}

		//write out? 
		if (output_filename!=NULL) {
			gzwrite(fi,buffer,read);
		}
	}
	
	//merge the results
	for (int i=0; i<threads; i++){ 
		ss->total_reads+=t_coverage[i];
		ss->skips+=t_skips[i];
		if (i>0) {
			for (int j=0; j<MAX_MEDIAN; j++) {
				t_median_counts[j]+=t_median_counts[i*MAX_MEDIAN+j];
			}
		} else {
			//cerr << genome_length << " " << total_entries << endl;
			if (total_entries<genome_length) {
				t_median_counts[0]=genome_length - total_entries;
			}
		}
		for (int j=0; j<bins; j++) {
			ss->gcbins[j]+=t_gcbins[i*bins+j];
		}
		for (int j=0; j<26; j++) {
			ss->reads_per_chr[j]+=t_reads_per_chr[i*26+j];
		}
	}

	//mean
	//ss->mean=((double)ss->total_reads)/entries;
	ss->mean=((double)ss->total_reads)/genome_length; //correct for full genome length

	//median
	int max_median_idx = 0;
	unsigned long reads_so_far=0;
	for (int j=0; j<MAX_MEDIAN; j++) {
		if (j<30) {
			cerr << "BIN DIST X " << j << " " << t_median_counts[j] << endl;
		}
		reads_so_far+=t_median_counts[j];
		if ( (total_entries/2)<reads_so_far ) {
			ss->median=j;
			break;
		}
	}

	//lets get a buffer to fit the file	
	//cout << "Done reading file  " << size_so_far <<  endl;

	free(t_reads_per_chr);
	free(t_skips);
	free(t_gcbins);
	free(t_median_counts);
	free(t_coverage);

	free(buffer);

	if (output_filename!=NULL) {
		gzclose(fi);	
	}

	return ss;
}

struct_stats * read_cov(char * filename, int bins, char ** ret_buffer, unsigned int * ret_entries, int * median) {
	size_t sz = 0 ; //size_so_far;	
	char * buffer = NULL;
	if (*ret_buffer == NULL ||  (filename!=NULL && strcmp("",filename)!=0)) {
		//cerr << "Reading coverage from file " << filename << endl;

		gzFile fptr = gzopen(filename,"r");
		if (fptr==NULL) {
			fprintf(stderr, "Failed to open file %s\n",filename);
			exit(1);
		}
		
		//find the size of one entry
		//size_t chunk = 30737418240L;
		size_t chunk = 1024*1024*1024L;
		size_t size_so_far = 0;
		buffer = (char*) malloc(chunk);
		if (buffer==NULL) {
			cerr << " FALLED TO MALLOC " << endl;
			exit(1);
		}


		//the real read loop
		while (!gzeof(fptr)) {
			size_t read = gzread(fptr,buffer+size_so_far,chunk);
			size_so_far+=read;
			buffer=(char*)realloc(buffer,size_so_far+chunk);
			if (buffer==NULL) {
				cerr << " FALLED TO REALLOC " << endl;
				exit(1);
			}
		}
		
		//lets get a buffer to fit the file	
		//cout << "Done reading file  " << size_so_far <<  endl;
		sz = size_so_far;
	} else {
		buffer = *ret_buffer;
		sz = *ret_entries*sizeof(struct_cov);
	}

	assert(sz%sizeof(struct_cov)==0);
	unsigned int entries = sz/sizeof(struct_cov);

	struct_stats * ss = find_stats(buffer, entries);	

	//cerr << "Entries " << entries << endl;
	*ret_buffer =buffer;
	*ret_entries = entries;

	return ss;
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
		genome_length+=lengths[i];
	}

	//compute the gcs
	gcs=(int**)malloc(sizeof(int*)*25);
	if (gcs==NULL) {
			cerr << "failed malloc for gcs" << endl;
			exit(1);
	}

	#pragma omp parallel for 
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
					ns++;
					break;
					//exit(1);
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
						ns++;
						break;
						//exit(1);
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
		free(fsta[i]);
	}
	free(buffer);
	free(ref);
	free(fsta);
	return NULL;	
}


void write_file(char * output_filename, char * buffer, unsigned int entries) {
	cerr << "writting file" << endl;
	gzFile *fi = (gzFile *)gzopen(output_filename,"wb");
        if (fi==NULL) {
                cerr<<"something terrible";
                return;
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


int main (int argc, char ** argv) {
	if (argc!=6) {
		cerr<<argv[0]<<" in_file_normal in_file_tumor ref output_tumor_filename output_normal_filename"<<endl;
		return 0;
	}

	omp_set_nested(1); 
	int threads=32;
	omp_set_num_threads(threads); //omp_get_num_threads();
	assert(sizeof(unsigned short)+sizeof(unsigned int)+sizeof(unsigned short)==sizeof(struct_cov));

	char * normal_coverage_filename=argv[1];
	char * tumor_coverage_filename=argv[2];
	char * fasta_filename=argv[3];
	char * output_tumor_filename = argv[4];
	char * output_normal_filename = argv[5];

	bins = 300;
	read_fasta(fasta_filename,bins);
	char * tumor_buffer = NULL;
	char * normal_buffer = NULL;
 	unsigned int tumor_entries, normal_entries;

	struct_stats * normal_stats = NULL;
	struct_stats * tumor_stats = NULL;

	cerr << "Reading in coverage files " << normal_coverage_filename << " and " << tumor_coverage_filename << endl;
	#pragma omp parallel num_threads(2)
	{
		int tidx=omp_get_thread_num();
		if (tidx==0) {
			normal_stats = stream_read_cov(normal_coverage_filename,bins, false, NULL, NULL, NULL, false);
		} else if (tidx==1) {
			tumor_stats = stream_read_cov(tumor_coverage_filename,bins, false, NULL, NULL, NULL, true);
		} else {
			cerr << " THREADING ERROR " << endl;
		}
	}
	print_stats(normal_stats);
	print_stats(tumor_stats);

	#pragma omp parallel num_threads(2)
	{
		int tidx=omp_get_thread_num();
		if (tidx==0) {
			//correct_coverage(normal_buffer, normal_entries, false, normal_stats, tumor_stats) ;
			//normal_stats = read_cov(NULL,bins, &normal_buffer, &normal_entries, &normal_median);
			//write_file(output_normal_filename,normal_buffer, normal_entries);
			normal_stats = stream_read_cov(normal_coverage_filename,bins, true,  normal_stats, tumor_stats, output_normal_filename, false);
		} else if (tidx==1) {
			//correct_coverage(tumor_buffer, tumor_entries, true, normal_stats, tumor_stats) ;
			//tumor_stats = read_cov(NULL,bins, &tumor_buffer, &tumor_entries, &tumor_median);
			//write_file(output_tumor_filename,tumor_buffer, tumor_entries);
			tumor_stats = stream_read_cov(tumor_coverage_filename,bins, true,  normal_stats, tumor_stats, output_tumor_filename, true);
		} else {
			cerr << " THREADING ERROR " << endl;
		}
	}
	print_stats(normal_stats);
	print_stats(tumor_stats);
	
}

//chr1    10000   N       3       AA^:A   D@+
