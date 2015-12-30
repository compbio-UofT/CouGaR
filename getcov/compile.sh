g++ -O3 get_cov.cpp -o get_cov -lz
g++ print_cov.cpp -o print_cov -lz -fopenmp -O3
g++ gc_normalize_cov.cpp -g -o normalize_cov -lz -fopenmp -O3 
g++ gc_normalize_stream_cov.cpp -g -o normalize_stream_cov -lz -fopenmp -O3 
g++ pickwindow.c igam.c -g -o pickwindow -O3
