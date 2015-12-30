g++ -Wall cnv_hmm.cpp  -lz -fopenmp -o hmm
g++ -Wall cnv_hmm_chrm.cpp  -lz -fopenmp -O3 -o hmm_chrm
g++ -Wall cnv_hmm_chrm.cpp  -lz -fopenmp -O3 -DNORMALIZED -o hmm_chrm_norm
