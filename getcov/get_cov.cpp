#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

using namespace std;

unsigned short get_chr(const char * s) {
	if (strlen(s)>3) {
		s+=3;
	}
	if (strlen(s)>0) {
		if (s[0]=='x' || s[0]=='X') {
			return 23;
		}
		if (s[0]=='y' || s[0]=='Y') {
			return 24;
		}
		if (s[0]=='m' || s[0]=='M') {
			return 25;
		}
		return atoi(s);
	} else {
		return 0;
	}
}

int main (int argc, char ** argv) {
	if (argc!=2) {
		cerr<<argv[0]<<" out_file";
		return 0;
	}

	char * filename=argv[1];

	/*FILE * fp =  fopen(filename,"wb");
	if (fp==NULL) {
		cerr<<"something terrible";
		return 1;
	}*/	

	char buffer[100];
	buffer[0]='\0';
	strcat(buffer, filename);
	strcat(buffer, ".gz");
	gzFile *fi = (gzFile *)gzopen(buffer,"wb");
	if (fi==NULL) {
		cerr<<"something terrible";
		return 1;
	}

	//cerr << buffer << endl;

	string chr,nucleo;
	unsigned int pos;
	unsigned short dp;
	string line;
	while(getline(cin, line)) {
		cin>>chr>>pos>>nucleo>>dp;
		//cout<<pos<<" "<<dp<<endl;
		unsigned short ichr=get_chr(chr.c_str());
		if (ichr==0 || ichr>25) {
			continue;
		} 
		gzwrite(fi,&ichr,sizeof(unsigned short));
		gzwrite(fi,&pos,sizeof(unsigned int));
		gzwrite(fi,&dp,sizeof(unsigned short));
	}

	gzclose(fi);


}

//chr1    10000   N       3       AA^:A   D@+
