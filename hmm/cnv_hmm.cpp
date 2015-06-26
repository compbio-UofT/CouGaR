#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <zlib.h>
#include <omp.h>
#include <queue>
#include <math.h>

#define THREADS	32



#define MIN_FLOW 6
#ifndef MAX_FLOW
#define MAX_FLOW 60
#endif
#define MAX_FREE 40
#define STATES 120


#define ZERO 1e-6


#define MAX_EDGE_SIZE 2400
//#define MAX_EDGE_SIZE 24000
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

using namespace std;

double tao = 0.0; //contamination

double abs(double a) {
	if (a<0) {
		return -a;
	}
	return a;
}

class pos {
	public:
		int chr;
		unsigned int coord;
		pos(int chr, unsigned int coord);
		pos();
		bool operator<(const pos &other) const;
		bool operator>(const pos &other) const;
		bool operator==(const pos &other) const;
		bool operator!=(const pos &other) const;
		string str() const;
};

class edge {
	public:
		pos posa;
		pos posb;
		edge(pos posa, pos posb);
		edge reverse();
		edge();
		unsigned int length() const;
		bool operator<(const edge &other) const;
		bool operator>(const edge &other) const;
		bool operator==(const edge &other) const;
		bool is_forward();
		unsigned int bound_cp();
};


class edge_info { 
	public:
		edge_info();
		unsigned int bp;
		double normal_coverage;
		double cancer_coverage;
		unsigned int type;
		int supporting;
		int copy_number;
};

//some global variables
set<pos> bps;

map<edge, edge_info > edges;
edge_info re_edges(edge & key) {
	if (edges.find(key)==edges.end()) {
		cerr << "ERROR IN EDGE LOOKUP" << endl;
		exit(1);
	}
	return edges[key];
}

map<pos, set<pos> > free_edges;
set<pos> re_free_edges(pos & key) {
	if (free_edges.find(key)==free_edges.end()) {
		cerr << "ERROR IN FREE EDGE LOOKUP" << endl;
		exit(1);
	}
	return free_edges[key];
}


unsigned long total_cancer_coverage=0;
unsigned long total_normal_coverage=0;
map<edge, int> free_edges_bound;
unsigned int total_normal_pair_arcs;
unsigned int total_cancer_pair_arcs;


unsigned int bp_range=100000;

edge fake_edge=edge(pos(0,0),pos(0,0));

//
// The pos class
//


pos::pos(int chr, unsigned int coord) {
	this->chr=chr;
	this->coord=coord;
}

pos::pos() {
	chr=-1;
	coord=-1;
}

bool pos::operator<(const pos &other) const {
	if (chr<other.chr) {
		return true;
	} else if (other.chr==chr) {
		return coord<other.coord;
	} else {
		return false;
	}
}

bool pos::operator>(const pos &other) const {
	if (chr>other.chr) {
		return true;
	} else if (other.chr==chr) {
		return coord>other.coord;
	} else {
		return false;
	}
}

bool pos::operator==(const pos &other) const {
	if (chr!=other.chr || coord!=other.coord) {
		return false;
	} 
	return true;

}


bool pos::operator!=(const pos &other) const {
	if (*this==other) {
		return false;
	}
	return true;

}

string pos::str() const {
	char buffer[5000];
	char cchr[10];
	if (chr==23) {
		cchr[0]='X'; cchr[1]='\0';
	} else if (chr==24) {
		cchr[0]='Y'; cchr[1]='\0';
	} else {
		sprintf(cchr,"%d",chr);
	}
	sprintf(buffer,"chr%s:%u" , cchr,coord);
	return string(buffer);
}


//
// The edge class
//
//if supporting > 0 then its a free edge!


bool edge::operator<(const edge &other) const {
	if (posa<other.posa) {
		return true;
	} else if (posa>other.posa) {
		return false;
	} else {
		return posb<other.posb;
	}
}

bool edge::operator>(const edge &other) const {
	if (posa>other.posa) {
		return true;
	} else if (posa<other.posa) {
		return false;
	} else {
		return posb>other.posb;
	}
}

edge edge::reverse() {
	return edge(posb,posa);
}

unsigned int edge::length() const {
	if (posa.chr==posb.chr) {
		if (posa.coord>posb.coord) {
			return posa.coord-posb.coord;
		} else {
			return posb.coord-posa.coord;
		}
	}
	return 0;
}

bool edge::operator==(const edge &other) const {
	return (posa==other.posa && posb==other.posb);
}

edge::edge(pos posa, pos posb) {
	this->posa=posa;
	this->posb=posb;
}

edge::edge() {
	*this=fake_edge;
}


bool edge::is_forward() {
	return posb>posa;
}



edge_info::edge_info() {
	bp=0;
	normal_coverage=0.0;
	cancer_coverage=0.0;
	type=5;
	supporting=0;
	copy_number=-1;
}





//take a chr name and give back the int
int to_chr(const char * s) {
	char buff[1024]="";
	unsigned int i=0;
	for (; i<strlen(s); i++) {
		buff[i]=tolower(s[i]);
	}
	buff[i]='\0';
	char * p = buff;
	if (i>3 && buff[0]=='c' && buff[1]=='h' && buff[2]=='r') {
		p=buff+3;
	}
	if (p[0]=='x') {
		return 23;
	}
	if (p[0]=='y') {
		return 24;
	}
	if (p[0]=='m') {
		return 25;
	}
	return atoi(p);
}

//read in the bp arcs
/*unsigned long read_arcs(char * filename, bool normal) {
	ifstream f (filename);
	string chr_s;
	unsigned int coord;
	int arcs;
	string line;
	getline(f,line);
	istringstream is(line);
	unsigned long total=0;
	is >> total;	
	if (total<1000) {
		cerr << "Warning probably read total wrong from arc links file!" << endl;
		exit(1);
	}

	cerr << "Total at head of file is " << total << endl;
	
	while (f) {
		getline(f,line);
		istringstream is(line);
		is >> chr_s >> coord >> arcs;
		int chr = to_chr(chr_s.c_str());
		
		pos p = pos(chr,coord);
		double v = ((double)arcs)/total;
		if (normal) {		 
			normal_pair_coverage[p]=v;
		} else {
			cancer_pair_coverage[p]=v;
		}
		//total+=arcs;
	}
	
	//set<pos>::iterator sit;
	//if (normal) {
	//	for (map<pos,double>::iterator sit = normal_pair_coverage.begin(); sit!=normal_pair_coverage.end(); sit++) {
	//		normal_pair_coverage[sit->first]=sit->second/total;
	//	}
	//	total_normal_pair_arcs=total;
	//} else {
	//	for (map<pos,double>::iterator sit = cancer_pair_coverage.begin(); sit!=cancer_pair_coverage.end(); sit++) {
	//		cancer_pair_coverage[sit->first]=sit->second/total;
	//	}
	//	total_cancer_pair_arcs=total;
	//}
	

	cerr << "Read " << total << " arcs from " << filename << endl;
	return total;
}*/

//read in the edges from clustering
void read_links(char * filename) {


	//insert the stoppers
	unsigned int lengths[]={249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,51304566,48129895,155270560,59373566,16571};
	for (int i=0; i<25; i++) {
		pos p = pos(i+1,MAX_EDGE_SIZE+1);
		bps.insert(p);
		p = pos(i+1,lengths[i]);
		bps.insert(p);
	}
		

	ifstream f (filename);
	//char chr[10]="";
	unsigned int bpa,bpb,l_from,l_to,cluster_idx;
	double avg_md;
	string chra,chrb;
	int type,total,total_links=0;
	while (f) {
		string line;
		getline(f,line);
		istringstream is(line);
		//chr2    144365053       185125297       3       2       2       0       2.8396e+08      1012    chr3    EDGE
		//chr3    84440059        190840923       3       9       67      70      3.04423e+08     1013    chr4    EDGE
		is >> chra >> bpa >> bpb >> type >> total >> l_from >> l_to >> avg_md >> cluster_idx >> chrb; 
		//cout << chra << bpa << bpb << type << total << l_from << l_to << avg_md << chrb << endl;
		//cout << line << endl << line << endl;
		int nchra=to_chr(chra.c_str());
		int nchrb=to_chr(chrb.c_str());
		if (bpa<=0 || bpb<=0) {
			continue;
		}
		total_links++;
		pos posa=pos(nchra,bpa);
		pos posb=pos(nchrb,bpb);

		if (bpa<MAX_EDGE_SIZE || bpb<MAX_EDGE_SIZE) {
			cerr << "SKIPPING LINK" << endl;
			continue;
		}

		if (nchra>26 || nchrb>26) {
			cerr << "ERROR READING IN LINKS" << endl;
			exit(1);
		}

		if (posa.chr>26 || posb.chr>26) {
			cerr << "EDGE EROR" << endl;
		}
		bps.insert(posa);
		bps.insert(posb);


		//double v = ((double)total)/total_paired;
		//if (v<1e-15) {
		//	cerr << "ERROR: READING LINKS " << v << endl;
		//	exit(1);
		//}
		

		
	
		//add the one direction
		edge ea = edge(posa,posb);

		free_edges[posa].insert(posb);
		edges[ea].type=type;
		edges[ea].supporting=total;
		edges[ea].bp=ea.length();

		//add the other direction
		edge eb = ea.reverse();

		if (type<=1) {
			type=1-type;
		}
		free_edges[posb].insert(posa);
		edges[eb].type=type;
		edges[eb].supporting=total;
		edges[eb].bp=eb.length();
	}
	cerr << "Read " << total_links << " links from " << filename << endl;	
	return ;
}




void read_cov(char * filename, bool normal) {
	cerr << "Reading coverage from file " << filename << endl;

	gzFile fptr = gzopen(filename,"r");
	if (fptr==NULL) {
		fprintf(stderr, "Failed to open file %s\n",filename);
		exit(1);
	}
	
	//find the size of one entry
	size_t soe = sizeof(unsigned short)+sizeof(unsigned int)+sizeof(unsigned short);
	size_t chunk = 30737418240L;
	//size_t chunk = 64*1024*1024;
	size_t size_so_far = 0;
	char * buffer = (char*) malloc(chunk);
	if (buffer==NULL) {
		cerr << " FALLED TO MALLOC " << endl;
		exit(1);
	}
	
	/*cerr << "WARNING" << endl;

	//the read shunt
	for (int i=0; i<2; i++) {
		size_t read = fread(buffer+size_so_far,1,chunk,fptr);
		size_so_far+=read;
		cerr << "Warning!!!!" << endl;
		buffer=(char*)realloc(buffer,size_so_far+chunk);
	}
	size_so_far=(size_so_far/soe)*soe;*/
	

	
	//the real read loop
	while (!gzeof(fptr)) {
		//size_t read = fread(buffer+size_so_far,1,chunk,fptr);
		size_t read = gzread(fptr, buffer+size_so_far,chunk); 
		size_so_far+=read;
		cerr << "Read so far " << size_so_far << endl;

		if (read==chunk) {
			cerr << "REALLOC" << endl;
			exit(1);
			buffer=(char*)realloc(buffer,size_so_far+chunk);
			if (buffer==NULL) {
				cerr << " FALLED TO REALLOC " << endl;
				exit(1);
			}
		}
	}
	
	//lets get a buffer to fit the file	
	cerr << "Done reading file  " << size_so_far <<  endl;

	size_t sz = size_so_far;	

	unsigned int entries = sz/soe;
	unsigned long total_coverage=0;
	unsigned long total_length=0;

	cerr << "started processing" << endl;

	
	//omp_set_num_threads(MIN(24,THREADS));
	#pragma omp parallel 
	{
	unsigned int threads = omp_get_num_threads();
	unsigned int thread_id = omp_get_thread_num();
	//cerr << "thread " <<  thread_id << endl;
	unsigned long total_coverage_t=0;
	unsigned long total_length_t=0;

	map<edge, edge_info > edges_t;

	set<pos>::iterator it = bps.begin();
	pos prev = *it; it++;
	//cout << prev.str() << endl;
	unsigned short chr, cov;
	unsigned int coord;
	//cerr << "started processing x2" << endl;
	for (unsigned int i=0; i<entries; i++) {
		if (i%threads!=thread_id) {
			continue; // not our job!
		}
		char* base = buffer+i*soe;
		chr=*((unsigned short *)base);
		if (chr==25 || chr==0) {
			continue;
		}
		base+=sizeof(unsigned short);
		coord=*((unsigned int *)base);
		base+=sizeof(unsigned int);
		cov=*((unsigned short *)base);

		total_coverage_t+=cov;

		
		pos p = pos(chr,coord);
		if (p.chr!=it->chr) {
			//cerr << "Switching chr from " << it->chr << " to " << p.chr << endl;
			//cerr << "Restarting looking for pos " << chr << " : " << coord << ", was at " << it->chr << " : " << it->coord << endl;
			//it=bps.lower_bound(p);
			//it--;
			it=bps.begin();
		}
		while (p>*it && it!=bps.end()) {
			prev=*it;
			it++;
		}
		//cout << p.str() << endl;
		if (it==bps.end()) {
			//cout << "broke at " << p.first << " " << p.second << "    " << (*it).first << " " << (*it).second << endl;
			break;	
		}
		if (prev.chr==chr) {
			edge ea = edge(prev,*it);
			edge eb = ea.reverse();
			total_length_t+=(*it).coord-prev.coord;

			if (normal) {
				edges_t[ea].normal_coverage+=cov;
				edges_t[eb].normal_coverage+=cov;
			} else {
				edges_t[ea].cancer_coverage+=cov;
				edges_t[eb].cancer_coverage+=cov;
			}
		}

	}

	#pragma omp critical 
	{
		//iterator over thread specific and add to main
		//cerr << " done processing critical " << thread_id << endl;
		total_coverage+=total_coverage_t;
		for (map<edge,edge_info>::iterator mit=edges_t.begin(); mit!=edges_t.end(); mit++) {
			edges[mit->first].normal_coverage+=mit->second.normal_coverage;
			edges[mit->first].cancer_coverage+=mit->second.cancer_coverage;
			total_length+=total_length_t;
		}
	}
	} //end openmp section
	cerr << " done processing " << endl;


	//lets take out the weird sections from the fractionization
	unsigned long total_coverage_effective = 0;
	double average = ((double)total_coverage)/((double)total_length);
	for (map<edge,edge_info>::iterator mit=edges.begin(); mit!=edges.end(); mit++) {

		if (normal && (mit->first).length()*average*5>=edges[mit->first].normal_coverage) {
			total_coverage_effective+=edges[mit->first].normal_coverage;
		}	
		
		if (!normal && (mit->first).length()*average*5>=edges[mit->first].cancer_coverage) {
			total_coverage_effective+=edges[mit->first].normal_coverage;
		}
	}



	for (map<edge,edge_info>::iterator mit=edges.begin(); mit!=edges.end(); mit++) {
	
		if (normal) {
			total_normal_coverage=total_coverage;
			edges[mit->first].normal_coverage/=total_coverage;
		} else {
			if (total_normal_coverage<100) {
				cerr << "Must read in normal before cancer coverages! " << endl;
				exit(1);
			}
			total_cancer_coverage=total_coverage;
			//edges[mit->first].cancer_coverage/=total_coverage; // old way without using tao
			double rt = edges[mit->first].cancer_coverage; //grab number of reads in tumor sample
			double rn = edges[mit->first].normal_coverage; // grab the fraction for normal, already a fraction!
			edges[mit->first].cancer_coverage=(rt/total_cancer_coverage - tao*rn)*(1/(1-tao));
		}
	}

	cerr << "total: " << total_coverage << endl;
	free(buffer);
}

int main(int argc, char ** argv) {
	//need to load in files
	if (argc!=5) {
		printf("%s links cov_cancer cov_normal contam(0-0.9)\n", argv[0]);
		exit(1);
	}



	char * links_filename=argv[1];
	char * cov_cancer_filename=argv[2];
	char * cov_normal_filename=argv[3];
	tao=atof(argv[4]);
	if (tao<0 || tao>0.9) {
		cerr << "contam " << tao << " is not in range (0-0.9) " << endl;
		exit(1);
	}

	cout << "#" << MAX_FLOW << "\t" << links_filename << "\t" << cov_cancer_filename << "\t" << cov_normal_filename << endl;

	//char * pairs_cancer_filename=argv[4];
	//char * pairs_normal_filename=argv[5];

	fake_edge=edge(pos(0,0),pos(0,0));

	//read in the arc weights
	//unsigned long total_normal_paired_mappings = read_arcs(pairs_normal_filename,true);	
	//unsigned long total_cancer_paired_mappings = read_arcs(pairs_cancer_filename,false);	
	//read in the free edges
	read_links(links_filename);

	cerr << "Slicing edges" << endl;
	//cerr << "Slicing WARNING edges" << endl;
	//lets slice up the rest
	set<pos> to_add;
	//unsigned int i=0;


	int last_chr=0;
	pos xprevious=pos(0,0);	
	for (set<pos>::iterator it=bps.begin(); it!=bps.end(); ) {
		pos current=*it;
		if (current.chr!=last_chr) {
			/*if (current.coord>MAX_EDGE_SIZE) {
				to_add.insert(pos(current.chr,current.coord-MAX_EDGE_SIZE));
			} else {
				cerr << "TOO CLOSE!" << current.str() << endl;
				exit(1);
			}*/
			to_add.insert(pos(current.chr,MAX(MAX_EDGE_SIZE,current.coord)-MAX_EDGE_SIZE));
			last_chr=current.chr;	
			if (xprevious.chr!=0) {
				to_add.insert(pos(xprevious.chr,xprevious.coord+MAX_EDGE_SIZE));
			}
		}
		/*i++;
		if (i>2000) {
			break;
		}*/
		//cout << current.chr << ":" << current.coord << endl;
		it++;
		pos next=*it;
		if (it==bps.end()) {
			break;
		}
		if (current.chr==next.chr && next.coord-current.coord>MAX_EDGE_SIZE) {
			int d = next.coord-current.coord;
			int e = d/(d%MAX_EDGE_SIZE==0 ? d/MAX_EDGE_SIZE : (d/MAX_EDGE_SIZE+1)); // d / (num requried edges )
			int k = e;
			if (d%e!=0) {
				d--;
				k++;
			}
			pos i = pos(current.chr,current.coord+k);
			while (next.coord>i.coord) {
				//cout <<  "added poss " << i.str() << "\t c " << current.str() << "\t n" << next.str() <<  endl;
				to_add.insert(i);
				k=e;
				if (d%e!=0) {
					d--;
					k++;
				}
				i = pos(current.chr,i.coord+k);
			}
		}
		xprevious=current;
	}
	if (xprevious.chr!=0) {
		to_add.insert(pos(xprevious.chr,xprevious.coord+MAX_EDGE_SIZE));
	}
	for (set<pos>::iterator it=to_add.begin(); it!=to_add.end(); it++) {
		if ( (*it).chr>26 ) {
			cerr << "ERROR IN SPLICING!" << endl;
			exit(1);
		}
		bps.insert(*it);
	}

	
	cerr << "Init edges" << endl;
	//initialize the rest of the edges, so can multithread
	for (set<pos>::iterator it=bps.begin(); it!=bps.end(); ) {
		free_edges[*it].size();		

		pos current=*it;
		//cout << current.chr << ":" << current.coord << endl;
		it++;
		pos next=*it;
		if (it==bps.end()) {
			break;
		}
	
		//else have current and next , lets figure it out
		if (current.chr==next.chr) {
			edge ea = edge(current,next);
			edges[ea].bp=ea.length();
			edge eb = ea.reverse();
			edges[eb].bp=eb.length();
		}
	}

	//read in the normal
	//cerr << "WARNING NO COV" << endl;
	read_cov(cov_normal_filename,true);
	read_cov(cov_cancer_filename,false);

	//read in the pairs coverage
	cerr << "Done reading in data..." << endl;

	set<pos>::iterator sit = bps.begin(); 
	pos previous;
	pos current=*sit;
	sit++;
	vector<edge> start_edges;
	while (sit!=bps.end()) {
		previous=current;
		current=*sit;
		if (current.chr==previous.chr) {
			edge ea = edge(previous,current);
			edge eb = ea.reverse();
			if (re_free_edges(ea.posa).size()>0 || re_free_edges(ea.posb).size()>0) {
				//cerr << " adding a start edge " << ea.posa.chr << ":" << ea.posa.coord << " " << ea.posb.chr << ":" << ea.posb.coord << endl;
				start_edges.push_back(ea);
				start_edges.push_back(eb);
			}
			//cerr << "Warning only forward edges " << endl;
		}
		sit++;
	}

	cerr << "Starting HMM..." << endl;

	double * states = (double*)malloc(sizeof(double)*STATES*(bps.size()+2));
	if (states==NULL) {
		cerr << "AMLLOC ERROC " << endl;
		exit(1);
	}

	
	for (int i=0; i<STATES; i++) {
		states[i]=1.0/STATES;
	}
	

	
	
	//lets just do chr1
	map<int,pair<edge,map<int,int> > > state_to_edge;
	int s=1;
	pos p = pos(0,0); 
	for (set<pos>::iterator it=bps.begin(); it!=bps.end(); it++) {
		pos c = *it;
		//cout << "C:" << c.str() << endl; 
		if (c.chr==p.chr) {
			edge e = edge(p,c);
			//cout << p.str() << "\t" << c.str() << endl; 
			edge_info ei = re_edges(e);
			unsigned long cancer_coverage = total_normal_coverage*ei.cancer_coverage/100;
			if (cancer_coverage==0) {
				cancer_coverage++;
			}
			unsigned long normal_coverage = (total_normal_coverage*ei.normal_coverage/100)/2;
			if (normal_coverage==0) {
				normal_coverage++;
			}
			//cout << normal_coverage << endl;		


	
			double emission[STATES];
			for (int i=0; i<STATES; i++) {
				if (normal_coverage>=30) {
					if (i==0) {
						emission[i]=-(((double)normal_coverage)*0.5)+cancer_coverage*log((((double)normal_coverage)*0.5));
					} else {
						emission[i]=-((double)normal_coverage*i)+cancer_coverage*log((normal_coverage*i));
					}
				} else {
					//TODO really should change this to just uniform...
					emission[i]=0;
				}
			}

			//if there is a free edge check if it lowers the copy count of increases
			set<pos> fs = re_free_edges(c);
			bool can_drop=false;	
			bool can_rise=false;
			for (set<pos>::iterator sit = fs.begin(); sit!=fs.end(); sit++) {
				edge e = edge(c,*sit);
				edge_info ei = re_edges(e);
				if (ei.type%2==0) {
					can_rise=true;
				} else {
					can_drop=true;
				}
			}
			double transistion[STATES*STATES];
			for (int i=0; i<STATES; i++) {
				double p=0;
				if (!can_rise && !can_drop) {
					double x = 0.999;
					for (int j=0; j<STATES; j++) {
						if (i==j) {
							p=log(x);
						} else {
							p=log((1-x)/(STATES-1));
						}
						transistion[STATES*i+j]=p;
					}
				}
				if (can_rise || can_drop) {
					double x = 0.6;
					for (int j=0; j<STATES; j++) {
						if (i==j) {
							p=log(x);
						} else {
							p=log((1-x)/(STATES-1));
						}
						transistion[STATES*i+j]=p;
					}
				}
				/*if (can_rise && !can_drop) {
					for (int j=0; j<5; j++) {
						if (i<=j) {
							p=log(0.99/(j-i+1));
						} else if (j>i) {
							p=log(0.0025/(5-(j-i+1)));
						}
						transistion[5*i+j]=p;
					}
				}
				if (!can_rise && can_drop) {
					for (int j=0; j<5; j++) {
						if (i>=j) {
							p=log(0.99/(i-j+1));
						} else if (j>i) {
							p=log(0.0025/(5-(i-j+1)));
						}
						transistion[5*i+j]=p;
					}
				}*/
			}



			map<int,int> back_t;

			for (int i=0; i<STATES; i++) {
				states[s*STATES+i]=transistion[STATES*i]+states[(s-1)*STATES+i]+emission[0];
				back_t[i]=0;
			}

			for (int i=0; i<STATES; i++) {
				for (int j=0; j<STATES; j++) {
					double  z = transistion[STATES*j+i]+states[(s-1)*STATES+j]+emission[i];
					if (z>states[s*STATES+i]) {
						states[s*STATES+i]=z;
						back_t[i]=j;
					}
				}	
			}		

			/*
			int min=0;
			for (int i=0; i<5; i++) {
				if (states[5*s+i]<states[5*s+min]) {
					min=i;
				}
			}
			for (int i=0; i<5; i++) {
				states[5*s+i]-=states[5*s+min];
			}*/

			state_to_edge[s]=pair<edge,map<int,int> >(e,back_t); 


			/*if (c.chr==1) {
				cout << p.str() << "\t" << c.str() << "\t" << ei.normal_coverage << "\t" << cancer_coverage << "\t" << normal_coverage << "\t" << ei.cancer_coverage/(ei.cancer_coverage*0.0001+ei.normal_coverage) << endl; 
				
				for (int i=0; i<5; i++) {
					cout << emission[i] << "," << states[s*5+i] << "\t";
				}
				int max=0;
				for (int i=0; i<5; i++) {
					if (states[s*5+i]>states[s*5+max]) {
						max=i;
					}
				}
				cout << "MAX " << max << endl;
				cout << endl;
			}*/
			s++;
		} else if (p.chr!=0) {
			//lets drop the states
			map<int, int> viterbi;
			int max=0; 
			s--;
			for (int i=0; i<STATES; i++) {
				if (states[s*STATES+i]>states[s*STATES+max]) {
					max=i;
				}
			}
			while (s>=0) {
				pair<edge, map<int,int> > & x = state_to_edge[s];
				edge & e = x.first;
				map<int,int> & back_t = x.second;
				//cout << max << "\t" << e.posa.str() << "\t" << e.posb.str() << endl;
				if (e.posa.chr!=0) {
					//cout << p.str() << "\t" << c.str() << endl; 
					re_edges(e);
					edges[e].copy_number=max;
					/*unsigned long cancer_coverage = total_normal_coverage*ei.cancer_coverage/100;
					unsigned long normal_coverage = (total_normal_coverage*ei.normal_coverage/100)/2;
					cout << e.posa.str() << "\t" << e.posb.str() << "\t" << ei.normal_coverage << "\t" << cancer_coverage << "\t" << normal_coverage << "\t" << ei.cancer_coverage/(ei.cancer_coverage*0.0001+ei.normal_coverage) << endl; */
				}
				max = back_t[max];
				s--;
			}
			for (int i=0; i<STATES; i++) {
				states[i]=1.0/STATES;
			}
			s=1;
			state_to_edge.clear();
				
		}
		p=c;
	}
		//lets drop the states
		map<int, int> viterbi;
		int max=0; 
		s--;
		for (int i=0; i<STATES; i++) {
			if (states[s*STATES+i]>states[s*STATES+max]) {
				max=i;
			}
		}
		while (s>=0) {
			pair<edge, map<int,int> > & x = state_to_edge[s];
			edge & e = x.first;
			map<int,int> & back_t = x.second;
			//cout << max << "\t" << e.posa.str() << "\t" << e.posb.str() << endl;
			if (e.posa.chr!=0) {
				//cout << p.str() << "\t" << c.str() << endl; 
				re_edges(e);
				edges[e].copy_number=max;
				/*unsigned long cancer_coverage = total_normal_coverage*ei.cancer_coverage/100;
				unsigned long normal_coverage = (total_normal_coverage*ei.normal_coverage/100)/2;
				cout << e.posa.str() << "\t" << e.posb.str() << "\t" << ei.normal_coverage << "\t" << cancer_coverage << "\t" << normal_coverage << "\t" << ei.cancer_coverage/(ei.cancer_coverage*0.0001+ei.normal_coverage) << endl; */
			}
			max = back_t[max];
			s--;
		}
		for (int i=0; i<STATES; i++) {
			states[i]=1.0/STATES;
		}
		s=1;
		state_to_edge.clear();
	

	edge e=fake_edge;
	int cp=-1;
	unsigned int normal=0;
	unsigned int cancer=0;
	for (set<pos>::iterator it = bps.begin(); it!=bps.end(); it++) {
		pos current = *it;
		if (e.posb.chr!=current.chr) {
			//starting a new chromosome
			if (e.length()>0) {
				//print the edge out
				cout << cp << "\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << e.length() << "\t" << normal << "\t" << cancer << endl;
			}
			e=edge(current,current);
			normal=0;
			cancer=0;
		} else {
			if (e.length()==0) {
				//initialize
				e=edge(e.posa,current);
				edge_info ei = re_edges(e);
				cancer=total_normal_coverage*ei.cancer_coverage/100;
				normal=(total_normal_coverage*ei.normal_coverage/100)/2;
				cp=re_edges(e).copy_number;
				if (re_free_edges(current).size()!=0) {
					cout << cp << "\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << e.length() << "\t" << cancer << "\t" << normal << endl;
					e=edge(current,current);
					normal=0;
					cancer=0;
				}
			} else {
				//see if we can add
				edge xe = edge(e.posb,current);
				int xcp = re_edges(xe).copy_number;
				if (xcp==cp) {
					//add it
					e.posb=current;
					edge_info xei = re_edges(xe);
					cancer+=total_normal_coverage*xei.cancer_coverage/100;
					normal+=(total_normal_coverage*xei.normal_coverage/100)/2;
					if (re_free_edges(current).size()!=0) {
						cout << cp << "\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << e.length() << "\t" << cancer << "\t" << normal << endl;
						e=edge(current,current);
						normal=0;
						cancer=0;
					}
				} else {
					//print and start again
					cout << cp << "\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << e.length() << "\t" << cancer << "\t" << normal << endl;
					e=edge(e.posb,current);
					edge_info ei = re_edges(e);
					cancer=total_normal_coverage*ei.cancer_coverage/100;
					normal=(total_normal_coverage*ei.normal_coverage/100)/2;
					cp=xcp;
					if (re_free_edges(current).size()!=0) {
						cout << cp << "\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << e.length() << "\t" << cancer << "\t" << normal << endl;
						e=edge(current,current);
						normal=0;
						cancer=0;
					}
				}
			}
		}
	}
	if (e.length()>0) {
		//print the edge out
		cout << cp << "\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << e.length() << "\t" << normal << "\t" << cancer << endl;
		//cout << cp << "\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << e.length() << endl;
	}

	cerr << "CLEAN" << endl;	
	return 0;
}
