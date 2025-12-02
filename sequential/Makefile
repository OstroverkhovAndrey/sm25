
all:
	mkdir -p build
	mkdir -p txt
	g++ -c prog.cpp -o  build/prog.o -I. -std=c++11 -O2 -fopenmp -g
	g++ -c variant.cpp -o build/variant.o -I. -std=c++11 -O2 -fopenmp -g
	g++ -c matrix.cpp -o build/matrix.o -I. -std=c++11 -O2 -fopenmp -g
	g++ -c args.cpp -o build/args.o -I. -std=c++11 -O2 -fopenmp -g
	g++ build/prog.o build/variant.o build/matrix.o build/args.o -o build/prog -std=c++11 -O2 -fopenmp -g
	./build/prog N 120 M 80 count_iter -1 delta 0.00001 --test