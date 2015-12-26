g++ -O3 get_cov.cpp -o get_cov -lz
g++ print_cov.cpp -o print_cov -lz -fopenmp
g++ gc_normalize_cov.cpp -o normalize_cov -lz -fopenmp
