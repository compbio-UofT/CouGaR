#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//\!f(k; \lambda)= \exp \left\{ {k\ln \lambda  - \lambda  - \ln \Gamma (k+1)} \right\},

double igam(double k, double m, double ep );

int main(int argc, char ** argv) {
	if (argc!=4) {
		fprintf(stdout, "%s read_arrival_rate genome_length cutoff\n" , argv[0]);
		exit(1);
	}
	const double arrival_rate = atof(argv[1]);
	const double genome_length = atof(argv[2]);
	const double cutoff = atof(argv[3]);


	double epsilon = 1.0;

	while ((1.0 + 0.5 * epsilon)!=1.0) {
	    epsilon = 0.5 * epsilon;
	}
	
	int window_size = 600;
	double ig=0.0;
	while (window_size<500000) {
		const double k = 1.5*arrival_rate*window_size; //arrival rate is assumed for diploid count
		const double m = window_size*arrival_rate;
		const double lnm = log(m);
		/*
		pdtrc=0;
		for (int kk=k; kk>=0; kk--) {
			pdtrc+=exp(kk*lnm-m-lgamma(kk+1));
		}
		pdtrc=1-pdtrc;*/
		ig = igam(k+1,m,epsilon)*genome_length;
		if (ig<cutoff) {
			break;
		}
		window_size+=100;
	}
	fprintf(stdout,"%e %d\n",ig,window_size);
	return 0;
}
