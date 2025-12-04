
all:
	mkdir -p build
	mkdir -p txt
	mpic++ -c prog.cpp -o  build/prog.o -I. -std=c++11 -O2 -fopenmp -g
	mpic++ -c variant.cpp -o build/variant.o -I. -std=c++11 -O2 -fopenmp -g
	mpic++ -c matrix.cpp -o build/matrix.o -I. -std=c++11 -O2 -fopenmp -g
	mpic++ -c args.cpp -o build/args.o -I. -std=c++11 -O2 -fopenmp -g
	mpic++ build/prog.o build/variant.o build/matrix.o build/args.o -o build/prog -std=c++11 -O2 -fopenmp -g
	mpirun -np 2 ./build/prog N 120 M 80 delta 0.00001 print_result