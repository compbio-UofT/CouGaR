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
#include <stdlib.h>
#include <string.h>
#include <limits>


//#define Q 1166516
#define Q 16

#define LARGE_COMPONENT 4

#define	SZ	400	//max ocpy count
#define MAX_WALK	180
#define MAX_REUSE	2
#define MIN_CP	1


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define	SOMATICW	1
#define GENOMICW	1

using namespace std;


int min_copies=0;
int multiplier=1;
int min_hmm_cp=2;


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
		bool genomic;
		edge(pos posa, pos posb,bool genomic);
		edge reverse();
		edge();
		unsigned int length() const;
		bool operator<(const edge &other) const;
		bool operator>(const edge &other) const;
		bool operator==(const edge &other) const;
		bool operator!=(const edge &other) const;
		bool is_forward();
		unsigned int bound_cp();
		edge canonical() const;
};


class edge_info { 
	public:
		edge_info();
		unsigned int length;
		unsigned int normal;
		unsigned int tumor;
		int type;
		int hmm_copy_number;
		int flow;
		double scores[SZ];
		double diffs[SZ];
		void poisson();
};


//some global variables
set<pos> bps;
map<pos,int> bp_support;

map<edge, edge_info > genomic_edges;
map<edge, edge_info > somatic_edges;
map<pos, set<pos> > jump_edges;

edge fake_edge=edge(pos(0,0),pos(0,0),true);

pos h1 = pos(0,0);
pos h2 = pos(0,1);
edge he = edge(h1,h2,true);

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
	} else if (chr==25) {
		cchr[0]='M'; cchr[1]='\0';
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
		if (posb==other.posb) {
			return genomic<other.genomic;
		} else {
			return posb<other.posb;
		}
	}
}

bool edge::operator>(const edge &other) const {
	if (posa>other.posa) {
		return true;
	} else if (posa<other.posa) {
		return false;
	} else {
		if (posb==other.posb) {
			return genomic>other.genomic;
		} else {
			return posb>other.posb;
		}
	}
}

edge edge::reverse() {
	return edge(posb,posa,genomic);
}

edge edge::canonical() const {
	if (posa<posb) {
		return edge(posa,posb,genomic);
	} else {
		return edge(posb,posa,genomic);
	}
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
	return (posa==other.posa && posb==other.posb && genomic==other.genomic);
}

bool edge::operator!=(const edge &other) const {
	return !(*this==other);
}

edge::edge(pos posa, pos posb,bool genomic) {
	this->posa=posa;
	this->posb=posb;
	this->genomic=genomic;
}

edge::edge() {
	*this=fake_edge;
}


bool edge::is_forward() {
	return posb>posa;
}



edge_info::edge_info() {
	length=0;
	normal=0;
	tumor=0;
	type=5;
	hmm_copy_number=-1;
	flow=0;
	for (int i=0; i<SZ; i++) {
		scores[i]=0;
		diffs[i]=0;
	}
}


void edge_info::poisson() {
		
	if (normal==0 || tumor==0) {
		//cerr << normal << " " << tumor << "NORMAL OR TUMOR IS ZERO, is this ok?" << endl;
	}


	// a bit of smoothing like thing?
	normal+=2;
	tumor+=2;

	//original formulation - CNVer
	for (int i=0; i<SZ; i++) {
		if (i==0 && min_copies==0) {
			scores[i]=normal*(min_copies==0 ? 0.1 : min_copies)-tumor*log(normal*(min_copies==0 ? 0.1 : min_copies));
		} else {
			scores[i]=normal*(min_copies+i*multiplier)-tumor*log(normal*(min_copies+i*multiplier));
		}
		double d0 = (min_copies==0 ? 0.1 : min_copies) +i*multiplier;
		double d1 = min_copies+(i+1)*multiplier;
		diffs[i]=normal*multiplier-tumor*log( d1/d0);
		diffs[i]=int(diffs[i]/10);
		//fprintf(stdout,"%d %lf\n",i,diffs[i]);
	}	
	//int doc_ratio=(tumor/normal)/multiplier;	
	double doc_ratio=(((double)tumor)/normal);	
	//the funky chicken formulation - not CNVer
	/*for (int i=doc_ratio; i<SZ; i++) {
		if (i==0 && min_copies==0) {
			scores[i]=normal*(min_copies==0 ? 0.1 : min_copies)-tumor*log(normal*(min_copies==0 ? 0.1 : min_copies));
		} else {
			scores[i]=normal*(min_copies+i*multiplier)-tumor*log(normal*(min_copies+i*multiplier));
		}
		double d0 = (min_copies==0 ? 0.1 : min_copies) +i*multiplier;
		double d1 = min_copies+(i+1)*multiplier;
		diffs[i]=normal*multiplier-tumor*log( d1/d0);
		diffs[i]=int(diffs[i]/10);
	}		
	for (int i=1; i<=doc_ratio; i++) {
		diffs[doc_ratio-i]=-diffs[doc_ratio+i];
	}
	if (doc_ratio==0 && tumor/normal>=1) {
		int real_doc_ratio=tumor/normal;
		double d0 = min_copies+real_doc_ratio;
		double d1 = min_copies+2*real_doc_ratio;
		diffs[0]=normal*real_doc_ratio-tumor*log( d1/d0);
		diffs[0]=-int(diffs[0]/10);
	}*/
	//end of the funky chicken
	/*for (int i=0; i<SZ; i++) {
		diffs[i]=(log(3+doc_ratio))*diffs[i];
	}*/
	//if (doc_ratio<=3.75) {
	if (doc_ratio<=5) {
		for (int i=0; i<SZ; i++) {
			diffs[i]=(length/2000)*multiplier;
		}
	}

	//make sure the min value is zero
	double min=10e100;
	for (int i=0; i<SZ; i++) {
		if (min>scores[i]) {
			min=scores[i];
		}
	}
	for (int i=0; i<SZ; i++) {
		//cerr << scores[i] << "\t" << scores[i]-min << endl;
		scores[i]-=min;
	}
	

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



map<edge, edge_info> edges;

void read_edges(char * filename) {
	FILE * fptr = fopen(filename,"r");
	if (fptr==NULL) {
		cerr << "failed to open file" << endl;
		exit(1);
	}
	char buffer[5000];
	double normal_coverage=0.0;
	double len=0.0;
	while (fgets(buffer,5000,fptr)) {
		if (buffer[0]=='#') {
			continue;
		}
		int cp,length,tumor,normal;
		long unsigned int coorda, coordb;
		char posa[100];
		char chra[10];
		int ichra=-1;
		char posb[100];
		char chrb[10];
		int ichrb=-1;
		int ret = sscanf(buffer,"%d\t%s\t%s\t%d\t%d\t%d",&cp,posa,posb,&length,&tumor,&normal);
		len+=length;
		normal_coverage+=normal;
	}
	fclose(fptr);

	double normal_average_coverage = normal_coverage/len;
	cerr << "# average normal coverage is " << normal_average_coverage << endl;


	fptr = fopen(filename,"r");
	if (fptr==NULL) {
		cerr << "failed to open file" << endl;
		exit(1);
	}
	
	while (fgets(buffer,5000,fptr)) {
		//fprintf(stdout,"#%s",buffer);
		//2       chrX:137812363  chrX:138585163  772800  264181  129674
		if (buffer[0]=='#') {
			continue;
		}
		int cp,length,tumor,normal;
		long int coorda, coordb;
		char posa[100];
		char chra[10];
		int ichra=-1;
		char posb[100];
		char chrb[10];
		int ichrb=-1;
		int ret = sscanf(buffer,"%d\t%s\t%s\t%d\t%d\t%d",&cp,posa,posb,&length,&tumor,&normal);
		int i=0;
		//cerr << buffer << endl;
		for (; posa[i]!=':'; i++) {
			chra[i]=posa[i];
		}
		chra[i]='\0';
		ichra=to_chr(chra);
		coorda=atol(posa+i+1);

		for (i=0; posb[i]!=':'; i++) {
			chrb[i]=posb[i];
		}
		chrb[i]='\0';
		coordb=atol(posb+i+1);
		ichrb=to_chr(chrb);

		//TODO
		int smoothing=10000;
			
		bool pass=false;
		double average_coverage = ((double)normal)/length;
		if (ichra==ichrb && ichra!=25 && average_coverage>10*normal_average_coverage) {
			cerr << "# Dropping region because cov " << average_coverage << " on " << chra << posa << " " << chrb << posb << " len " <<  length << " , normal avg cov is " << normal_average_coverage << endl;
			normal=0;
			tumor=0;
			bool pass=true;
		}

		//TRY TO FIX THE HMM PROBLEM : HMM ASSUMED DIPLOID CHRM .. so normal should be expected per single copy, multiplying by 2 here should undo this.. i know ugly.... but if it works... then that would be great... 
		if (ichra==ichrb && ichra==25) {
			normal*=2;
		}



		if (cp>min_hmm_cp || (ichra==ichrb && ( coordb-coorda<smoothing || pass))) {
			//insert the nodes
			if (ichra==ichrb && coordb-coorda>2) {
				pos from = pos(ichra,coorda);
				pos middle = pos(ichra,(coorda+coordb)/2);
				pos to = pos(ichrb,coordb);

				bps.insert(from);
				bps.insert(middle);
				bps.insert(to);


				//add the edges
				edge efm = edge(from,middle,true);
				edge emt = edge(middle,to,true);
		
				edge_info ei;
				ei.length=length/2;
				ei.type=0;
				ei.hmm_copy_number=cp;

				ei.normal=normal/2;
				ei.tumor=tumor/2;
				ei.poisson();	
				//190231633	
				/*if (coorda==121356294 ) {
					cerr << "XXX" << posa << " " << posb << " " << ei.tumor << " " << ei.normal << " "  << ei.tumor/ei.normal << " " << ei.diffs[0] << " " << ei.diffs[1] << endl;
					exit(1);
				}*/
		


				genomic_edges[efm]=ei;
				genomic_edges[efm.reverse()]=ei;
				genomic_edges[emt]=ei;
				genomic_edges[emt.reverse()]=ei;
			} else {
				pos from = pos(ichra,coorda);
				pos to = pos(ichrb,coordb);

				bps.insert(from);
				bps.insert(to);


				//add the edges
				edge e = edge(from,to,true);
		
				edge_info ei;
				ei.length=length;
				ei.type=0;
				ei.hmm_copy_number=cp;

				ei.normal=normal;
				ei.tumor=tumor;
				ei.poisson();

				genomic_edges[e]=ei;
				genomic_edges[e.reverse()]=ei;

			}
			
		}
		
	}
}


void read_links(char * filename) {
	FILE * fptr = fopen(filename,"r");
	if (fptr==NULL) {
		cerr << "failed to open file" << endl;
		exit(1);
	}
	char buffer[5000];
	while (fgets(buffer,5000,fptr)) {
		if (buffer[0]=='#') {
			continue;
		}
		//chr1    59446789        +       chr1    59446926        +       13
		char chra[10];
		int ichra=-1;
		char stranda;
		unsigned long coorda;

		char chrb[10];
		int ichrb=-1;
		char strandb;
		unsigned long coordb;

		int support;


		int ret = sscanf(buffer,"%s\t%lu\t%c\t%s\t%lu\t%c\t%d",chra,&coorda,&stranda,chrb,&coordb,&strandb,&support);

		edge_info ei;

		ichra=to_chr(chra);
		ichrb=to_chr(chrb);

		if (ichra>25 || ichrb>25 || ichra==0 || ichrb==0) {
			cerr << " skipping link, chrM or chr?? " << endl;
			continue;
		}

		pos posa=pos(ichra,coorda);
		if (bps.find(posa)==bps.end()) {
			continue;
		}	
		pos posb=pos(ichrb,coordb);
		if (bps.find(posb)==bps.end()) {
			continue;
		}	

		//get the edge type
		int type=0;
		if (stranda=='+') {
			if (strandb=='+') {
				type=0;	
			} else {
				//negative
				type=2;
			}
		} else {
			if (strandb=='+') {
				type=3;
			} else {
				type=1;
			}
		}


		//get the normal support
		if (bp_support.find(posa)==bp_support.end() ) {
			//cerr << "failed to find posa in bp_support " << endl;
		}
		if (bp_support.find(posb)==bp_support.end() ) {
			//cerr << "failed to find posb in bp_support " << endl;
		}
		int normal=(bp_support[posa]+bp_support[posb])/4; //the lambda , average of two sites (2 copies each)
		int tumor=support;	

		edge e = edge(posa,posb,false);
		if (somatic_edges.find(e)!=somatic_edges.end() && somatic_edges[e].tumor>=ei.tumor) {
			//skip this already have a better one
			//cerr << "SKipping edge " << buffer << endl;
			continue;
		}

		//set up the ei
		ei.length=0;
		ei.type=type;
		//cerr << posa.str() << " " << posb.str() << " " << normal << " " << tumor << endl;	
		ei.normal=normal;
		ei.tumor=tumor;
		ei.poisson();
		

		//put in the forward
		somatic_edges[e]=ei;

		//put in the reverse
		if (type<2) {
			type=1-type;
		}
		somatic_edges[e.reverse()]=ei;
		somatic_edges[e.reverse()].type=type;

		

		jump_edges[posa].insert(posb);
		jump_edges[posb].insert(posa);

		
	}
}




void read_bp_coverages(char * filename) {
	FILE * fptr = fopen(filename,"r");
	if (fptr==NULL) {
		cerr << "failed to open file" << endl;
		exit(1);
	}
	char buffer[5000];
	while (fgets(buffer,5000,fptr)) {
		if (buffer[0]=='#') {
			continue;
		}
		//chr1    59446926        12
		char chra[10];
		int ichra=-1;
		unsigned long coorda=0;
		int support=0;
		int ref = sscanf(buffer,"%s\t%lu\t%d",chra,&coorda,&support);
		ichra=to_chr(chra);
		bp_support[pos(ichra,coorda)]=support;
	}

}


set<pos> get_genomic_edges(pos p) {
	set<pos> r;
	set<pos>::iterator sit=bps.find(p);
	if (sit==bps.end()) {
		cerr << " A BAD ERROR IN GET GENOMIC EDGES " << endl;
		exit(1);
	}	

	//try the previous position
	sit--;
	if (sit!=bps.end() ) {
		pos previous = *sit;
		edge e = edge(previous, p,true);
		if (genomic_edges.find(e)!=genomic_edges.end()) {
			r.insert(previous);
		}
	}
	
	//try the next position
	sit=bps.find(p);
	sit++;
	if (sit!=bps.end() ) {
		pos next = *sit;
		edge e = edge(p, next,true);
		if (genomic_edges.find(e)!=genomic_edges.end()) {
			r.insert(next);
		}
	}

	return r;

	
}


void find_connected_component(pos p, set<pos> & so_far) {
	//add this node to the set
	if (so_far.count(p)>0) {
		cerr << " An error has happened" << endl;
		exit(1);
	}
	so_far.insert(p);

	//try going out on all genomic edges
	set<pos> ges = get_genomic_edges(p);
	for (set<pos>::iterator xit = ges.begin(); xit != ges.end(); xit++) {
		pos x = *xit;
		if (so_far.count(x)==0) {
			find_connected_component(x,so_far);
		}
	}

	//try going out on all free edges	
	set<pos> jumps = jump_edges[p];
	for (set<pos>::iterator xit = jumps.begin(); xit != jumps.end(); xit++) {
		pos x = *xit;
		if (so_far.count(x)==0) {
			find_connected_component(x,so_far);
		}
	}		
}


map<pos, int> find_connected_components() {
	int cid=0;
	map<pos, int> connected;
	for (set<pos>::iterator sit = bps.begin(); sit!=bps.end(); sit++) {
		pos p = *sit;
		//already found this component
		if (connected.find(p)!=connected.end()) {
			continue;
		}

		//find the the component
		set<pos> component;
		find_connected_component(p,component);
		for (set<pos>::iterator xit = component.begin(); xit != component.end(); xit++) {
			pos x = *xit;
			if (connected.find(x)!=connected.end()) {
				cerr << "double found a component!" << endl;
				exit(1);
			}
			connected[x]=cid;
		}
		/*if (component.size()>2) {	
			cout << "found component of size " << component.size() << endl;
		}*/
		cid++;
	}
	
	return connected;

}





string edges_summary(map<edge,edge_info> & m) {
	stringstream ss;
	for (map<edge,edge_info>::iterator mit = m.begin(); mit!=m.end(); mit++) {
		const edge & e = mit->first;
		if (e!=e.canonical()) {
			continue;
		}
		const edge_info & ei = mit->second;
		ss << "(" << e.posa.str() << "," << e.posb.str() << " , hmm_cp: " << ei.hmm_copy_number << " flow: " << ei.flow  << " , cancer: " << ei.tumor << " , normal: " << ei.normal <<  " , ratio: " << ei.tumor/(1+ei.normal)  << endl;
	}
	return ss.str();
}

string arc_strings(int from_node, int to_node , int type, int low, int cap, int cost) {
	stringstream ss;

	char stranda='+';
	char strandb='+';
	if (type==0 || type==2) {
		stranda='+';
		if (type==0) {
			strandb='+';
		} else {
			strandb='-';
		}
	} else {
		stranda='-';
		if (type==1) {
			strandb='-';
		} else {
			strandb='+';
		}
	}


	//if ( (to_node==2 && from_node!=1) || (from_node==3 && to_node!=1) || (from_node==1 && to_node!=2)) {
	//	ss << "a\t" << (stranda=='+' ? -1 : 0 )+2*from_node << "\t" << (strandb=='+' ? -1 : 0 )+2*to_node << "\t" << 0 << "\t" << 0 << "\t" << cost << endl;  
	//} else {
		ss << "a\t" << (stranda=='+' ? -1 : 0 )+2*from_node << "\t" << (strandb=='+' ? -1 : 0 )+2*to_node << "\t" << (to_node==2 ? low : 0) << "\t" << cap << "\t" << cost << endl; 
	//} 

	//if ((from_node==2 && to_node!=1) || (to_node==3 && from_node!=1) || (to_node==1 && from_node!=2)) {
	//	ss << "a\t" << (strandb=='-' ? -1 : 0 )+2*to_node << "\t" << (stranda=='-' ? -1 : 0 )+2*from_node << "\t" << 0 << "\t" << 0 << "\t" << cost << endl;  
	//} else {
		ss << "a\t" << (strandb=='-' ? -1 : 0 )+2*to_node << "\t" << (stranda=='-' ? -1 : 0 )+2*from_node << "\t" << (from_node==2 ? low : 0) << "\t" << cap << "\t" << cost << endl;  
	//}
	return ss.str();
}

void add_flow_to_edge(pos posa,pos posb) {
	//lets find the cheapeast way for this flow
	double somatic_cost=10e100;
	double genomic_cost=10e100;
	
	bool has_genomic=false;
	edge ge = edge(posa,posb,true);
	if (genomic_edges.find(ge)!=genomic_edges.end()) {
		has_genomic=true;
	}

	bool has_somatic=false;
	edge se = edge(posa,posb,false);
	if (somatic_edges.find(se)!=somatic_edges.end()) {
		has_somatic=true;
	}

	if (has_genomic) {
		edge_info gei=genomic_edges[ge];
		if ((gei.flow/2)<SZ-1) {
			genomic_cost=gei.scores[gei.flow/2+1]-gei.scores[gei.flow/2];
			//cout << "gflow is " << gei.flow << endl;
			//cerr << "GEN COST " << genomic_cost << endl;
		}
	}

	if(has_somatic) {
		edge_info sei = somatic_edges[se];
		if ((sei.flow/2)<SZ-1) {
			somatic_cost=sei.scores[sei.flow/2+1]-sei.scores[sei.flow/2];
			//cout << "sflow is " << sei.flow << endl;
			//cerr << "SOM COST " << somatic_cost << endl;
		}
	}
	
	if (somatic_cost>10e99 && genomic_cost>10e99) {
		//cerr << "SOMTHIGN BAD" << endl;	
		cerr << "OVERFLOW!\t" << posa.str() << "\t" << posb.str() << "\t" << (has_genomic ? genomic_edges[ge].flow : -1 ) << "\t" << (has_somatic ? somatic_edges[se].flow : -1 ) << endl;
		
		//exit(1);
	} else {
		if (somatic_cost<genomic_cost) {
			somatic_edges[se].flow++;
			somatic_edges[se.reverse()].flow++;
		} else {
			genomic_edges[ge].flow++;
			genomic_edges[ge.reverse()].flow++;
		}
	}

}

void flow_solve() {
	cerr << "Running flow with ? contigs" << endl;
	//output the graph!
	stringstream ss;
	int num_nodes=0;
	map<pos,int> node_ids;
	map<int,pos> inv_node_ids;

	num_nodes++;
	pos source = pos(0,1);
	node_ids[source]=num_nodes;
	inv_node_ids[num_nodes]=source;
	
	num_nodes++;
	pos sourcex = pos(0,10);
	node_ids[sourcex]=num_nodes;
	inv_node_ids[num_nodes]=sourcex;

	num_nodes++;
	pos sink = pos(0,2);
	node_ids[sink]=num_nodes;
	inv_node_ids[num_nodes]=sink;

	for (set<pos>::iterator sit=bps.begin(); sit!=bps.end(); sit++ ) {
		num_nodes++;
		node_ids[*sit]=num_nodes;
		inv_node_ids[num_nodes]=*sit;
	}
	num_nodes*=2;
	for (int i=1; i<=num_nodes; i++) {
		ss << "c NODE " << i << " " << inv_node_ids[(i+1)/2].str() << endl;
		ss << "n\t" << i << "\t0" << endl; 
	}


	//print the s and t edges
	int arcs=0;
	int low=0;
	int cap=0;
	int cost=0;
	
	ss << arc_strings(1,2,0,0,SZ,Q);
	arcs+=2;
	/*ss << "a\t" << 1 << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\t" << 100 << endl;  
	arcs+=1;
	ss << "a\t" << 2 << "\t" << 4 << "\t" << 0 << "\t" << 0 << "\t" << Q << endl;  
	arcs+=1;*/

	for (int i=4; i<(num_nodes/2); i++) {
		int low=0;
		//int cap=contigs;
		int cost=0;
		int cap=0;
		ss << arc_strings(2,i,0,0,SZ,0);
		ss << arc_strings(2,i,2,0,SZ,0);
		arcs+=4;
		ss << arc_strings(i,1,0,0,SZ,0);
		ss << arc_strings(i,1,3,0,SZ,0);
		arcs+=4;
		/*ss << "a\t" << 4 << "\t" << 2*i << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;  
		ss << "a\t" << 4 << "\t" << 2*i-1 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;  
		arcs+=2;
		ss << "a\t" << 2*i << "\t" << 2 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;  
		ss << "a\t" << 2*i-1 << "\t" << 2 << "\t" << 0 << "\t" << 0 << "\t" << 0 << endl;  
		arcs+=2;*/
	}
	
	ss << "c Here are the genomic edges" << endl;
	for (map<edge,edge_info>::iterator mit = genomic_edges.begin(); mit!=genomic_edges.end(); mit++) {
		const edge & e = mit->first;
		if (e!=e.canonical()) {
			continue;
		}
		mit->second.flow=0;
		edge_info & ei = mit->second;
		if (node_ids.find(e.posa)==node_ids.end()) {
			cerr << "Failed to find something ... " << endl;
			exit(1);
		}
		if (node_ids.find(e.posb)==node_ids.end()) {
			cerr << "Failed to find something ... " << endl;
			exit(1);
		}
		for (int i=0; i<SZ; ) {
			int low=0;
			int cap=0;
			int cost=0;
			if (i>0) {
				//cost=ei.scores[i]-ei.scores[i-1];
				cost=ei.diffs[i-1];
				while (i<SZ  && (cost==(int)(ei.diffs[i-1]))) {
					i++;
					cap++;
				}			
			} else {
				i++;
			}

			if (cap!=0) {
				if (cost>-2000000000L) {
					ss << "c Genomic\t" << e.posa.str() << "\t" << e.posb.str()  << "\t" << cost << "\t" << cap << "\t" << e.length() << "\t" << ei.normal << "\t" << ei.tumor << endl;
					ss << arc_strings(node_ids[e.posa],node_ids[e.posb],0,low,cap,cost);
					arcs+=2;
				}
			}

		}
	}

	ss << "c Here are the somatic edges" << endl;	
	for (map<edge,edge_info>::iterator mit = somatic_edges.begin(); mit!=somatic_edges.end(); mit++) {
		const edge & e = mit->first;
		if (e!=e.canonical()) {
			continue;
		}
		mit->second.flow=0;
		edge_info & ei = mit->second;
		if (node_ids.find(e.posa)==node_ids.end()) {
			cerr << "Failed to find something ... " << endl;
			exit(1);
		}
		if (node_ids.find(e.posb)==node_ids.end()) {
			cerr << "Failed to find something ... " << endl;
			exit(1);
		}
		for (int i=0; i<SZ; ) {
			int low=0;
			int cap=0;
			int cost=0;
			if (i>0) {
				//cost=ei.scores[i]-ei.scores[i-1];
				cost=ei.diffs[i-1];
				while (i<SZ  && (cost==(int)(ei.diffs[i-1]))) {
					i++;
					cap++;
				}			
			} else {
				i++;
			}

			//cost+=1000*multiplier; //TODO

			//cost=10*multiplier; //TODO HARD CODED	
			cost=((int)Q/0.9)*multiplier;
			if (cap!=0) {
				ss << "c Somatic\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << ei.type << "\t" << cost << "\t" << cap << "\t" << ei.normal << "\t" << ei.tumor << endl;
				ss << arc_strings(node_ids[e.posa],node_ids[e.posb],ei.type,low,cap,cost);
				arcs+=2;
			}
		}
	}
	


	cout << "c Here goes nothing ... " << endl;
	cout << "p\tmin\t" << num_nodes << "\t" << arcs << endl;
	cout << ss.str() << endl;

	/*
	//open the problem file
	gzFile gzout = gzopen("./problem_file.gz", "wb");
	char *b;
	b=(char*)malloc(sizeof(char)*(1000+ss.str().size()));
	if (b==NULL) {
		fprintf(stderr, "Failed to malloc for problem file output \n");
		exit(1);
	}
	
	sprintf(b,"c Here goes nothing ... \np\tmin\t%d\t%d\n%s" , num_nodes,arcs, ss.str().c_str());
	gzwrite(gzout,b,strlen(b)); 
	gzclose(gzout);
	*/
	//ofstream fs ("./problem_file");
	//fs << "c Here goes nothing ... " << endl;
	//fs << "p\tmin\t" << num_nodes << "\t" << arcs << endl;
	//fs << ss.str();
	//fs.close();

	//run cs2.exe on it
	/*
	FILE *fp;
	int status;
	
	fp = popen("zcat ./problem_file.gz | grep -v \"^c\" | /data/misko/2013.04.12/cs2-4.3/cs2.exe", "r");
	if (fp == NULL) {
		printf("Failed to run cs2.exe command\n" );
		exit(1);
	}

	char buffer[1035];
	while (fgets(buffer, sizeof(buffer)-1, fp) != NULL) {
		//printf("SOLUTION %s", buffer);
		switch(buffer[0]) {
			case 'c':
				cout << buffer;
				break;
			case 's':
				break;
			case 'f':
				int from,to,ret,f,fromb,tob;
	 			ret = sscanf(buffer,"f\t%d\t%d\t%d\n",&from,&to,&f);
				if (ret!=3) {
					cerr << "Failed to parse flow output " << endl;
					exit(1);
				}
				fromb=(1+from)/2;
				tob=(1+to)/2;
				if ( (fromb>3 && inv_node_ids.find(fromb)==inv_node_ids.end()) || (tob>3 && inv_node_ids.find(tob)==inv_node_ids.end() )) {
					cerr << "Failed to inverse lookup flow nodes " << endl;
					exit(1);
				}

				//it actual genomic or somatic edge
				if (f>0 && fromb>3 && tob>3) {
					//if ( (tob==58 && fromb==47) || (tob==47 && fromb==58) ) {
						//printf("%s", buffer);
						//if (f>1) {
						//	cerr << "only expecting unit flow!" << endl;
						//	exit(1);
						//}
						pos posa = inv_node_ids[fromb];
						pos posb = inv_node_ids[tob];
						//cout << posa.str() << "\t" << from << "\t" << posb.str() << "\t" << to << endl;
						while (f>0) {
							add_flow_to_edge(posa,posb);
							f--;
						}
					//}
				}
				break;
			default:
				cerr << "failed o hhandl this " << endl;
				exit(1);
		}
	}

	//correct the flow counts
	for (map<edge,edge_info>::iterator mit=genomic_edges.begin(); mit!=genomic_edges.end(); mit++) {
		mit->second.flow/=2;
	}
	for (map<edge,edge_info>::iterator mit=somatic_edges.begin(); mit!=somatic_edges.end(); mit++) {
		mit->second.flow/=2;
	}
	
	cout << "GENOMIC" << endl;
	cout << edges_summary(genomic_edges);
	
	cout << "SOMATIC" << endl;
	cout << edges_summary(somatic_edges);

	pclose(fp);*/

}

int main ( int argc, char ** argv) {
	if (argc!=7) {
		cerr << argv[0] << " links edges [only large component?Y/N [N]] min_copies[0] multiplier[1] min_hmm_cp[2]" << endl;
		exit(1);
	}

	
	//output the command line
	cerr << "#CMD-LINE: " ;
	for (int i=0; i<argc-1; i++) {
		cerr << argv[i] << " ";
	}
	cerr << argv[argc-1] << endl;

	bool only_largest=false;	
	if (argv[3][0]=='Y' || argv[3][0]=='y') {
		cerr << "# only largest component" << endl;
		only_largest=true;
	} else {
		cerr << "# all components" << endl; 
		only_largest=false;
	}

	min_copies=atoi(argv[4]);
	if (min_copies<0 || min_copies>1000) {
		cerr << "ERROR in min copies (range 0-1000) " << endl;
		exit(1);
	}

	multiplier=atoi(argv[5]);
	if (multiplier<1 || multiplier>1000) {
		cerr << "ERROR in multiplier (range 1-1000) " << endl;
		exit(1);
	}

	min_hmm_cp=atoi(argv[6]);
	if (min_hmm_cp<0 || min_hmm_cp>2000) {
		cerr << "ERROR in min_hmm_cp (range 0,2000) \n";
		exit(1);
	}

	char * links_filename = argv[1];
	//char * bp_coverages_filename = argv[2];
	char * edges_filename = argv[2];


	read_edges(edges_filename);
	//read_bp_coverages(bp_coverages_filename);
	read_links(links_filename);
	fprintf(stderr, "# DONE READ IN!\n");
	/*map<pos,int> connected_components = find_connected_components();
	map<int, int> connected_components_sizes;
	int size_3_or_more=0;
	set<int> large_components;
	for (map<pos,int>::iterator mit = connected_components.begin(); mit!=connected_components.end(); mit++) {
		connected_components_sizes[mit->second]++;
	}

	//lets find the max comopnent and look at that
	int max_id=-1;
	int max=-1;
	map<int, int> sizes;
	for (map<int, int>::iterator mit = connected_components_sizes.begin(); mit!=connected_components_sizes.end(); mit++) {
		if (mit->second>2) {
			size_3_or_more++;
		}
		if (mit->second>LARGE_COMPONENT) {
			large_components.insert(mit->first);
		}
		sizes[mit->second]+=1;
		if (mit->second>max) {
			max_id=mit->first;
			max=mit->second;
		}
	}	
		
	if (max_id==-1) {
		cerr << "No max components found....";
		exit(1);
	}

	cout << "c SOMATICW: " << SOMATICW << "\t" << "GENOMICW: " << GENOMICW << endl;
	cerr << "Found " << connected_components_sizes.size() << " connected components, >2 " << size_3_or_more <<  endl;
	for (map<int,int>::iterator mit = sizes.begin(); mit!=sizes.end(); mit++) {
		cerr << "SZ: " << mit->first << " " << mit->second << endl;
	}
	//exit(1);
	cerr << "Considering max component of size " << max << " id " << max_id << endl;

	//drop everything but this component
	//genomic first
	set<edge> to_remove;
	for (map<edge, edge_info>::iterator mit = genomic_edges.begin(); mit!=genomic_edges.end(); mit++ ) {
		edge e = mit->first;
		if (only_largest) {
			if (connected_components[e.posa]!=max_id || connected_components[e.posb]!=max_id) {
				to_remove.insert(e);
			} 
		} else {
			if (large_components.count(connected_components[e.posa])==0 || large_components.count(connected_components[e.posb])==0) {
				to_remove.insert(e);
			}
		}
	}
	for (set<edge>::iterator sit = to_remove.begin(); sit!=to_remove.end(); sit++ ){
		genomic_edges.erase(*sit);
	}
	to_remove.clear();
	//now drop the somatic
	for (map<edge, edge_info>::iterator mit = somatic_edges.begin(); mit!=somatic_edges.end(); mit++ ) {
		edge e = mit->first;
		if (only_largest) {
			if (connected_components[e.posa]!=max_id || connected_components[e.posb]!=max_id) {
				to_remove.insert(e);
			}
		} else {
			if (large_components.count(connected_components[e.posa])==0 || large_components.count(connected_components[e.posb])==0) {
				to_remove.insert(e);
			}
		}
	}
	for (set<edge>::iterator sit = to_remove.begin(); sit!=to_remove.end(); sit++ ){
		somatic_edges.erase(*sit);
	}
	//reconstruct bps
	bps.clear();
	for (map<pos,int>::iterator mit=connected_components.begin(); mit!=connected_components.end(); mit++ ){
		if (only_largest) {
			if (mit->second==max_id) {
				bps.insert(mit->first);
			}
		} else {
			if (large_components.count(mit->second)!=0) {
				bps.insert(mit->first);
			}
		}
	}*/

	//now we just have the largest component left!
	cerr << "SOLVING FLOW FOR " << endl;
	flow_solve();
	
	return 0;
}
