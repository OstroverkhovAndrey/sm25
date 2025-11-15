
all:
	mkdir -p build
	mpic++ -c prog.cpp -o  build/prog.o -I. -std=c++11 -O2 -fopenmp
	mpic++ -c variant.cpp -o build/variant.o -I. -std=c++11 -O2 -fopenmp
	mpic++ -c matrix.cpp -o build/matrix.o -I. -std=c++11 -O2 -fopenmp
	mpic++ -c args.cpp -o build/args.o -I. -std=c++11 -O2 -fopenmp
	mpic++ build/prog.o build/variant.o build/matrix.o build/args.o -o build/prog -std=c++11 -O2 -fopenmp
	mpirun -np 2 ./build/prog N 60 M 40 count_iter 10 delta 0.00001 --test