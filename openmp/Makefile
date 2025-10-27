
all:
	mkdir -p build
	g++ -c prog.c -o  build/prog.o -I. -std=c++20 -O2 -fopenmp
	g++ -c variant.c -o build/variant.o -I. -std=c++20 -O2 -fopenmp
	g++ -c matrix.c -o build/matrix.o -I.  -O2 -fopenmp
	g++ -c args.c -o build/args.o -I. -O2 -fopenmp
	g++ build/prog.o build/variant.o build/matrix.o build/args.o -o build/prog -O2 -fopenmp
	./build/prog N 1200 M 800 count_iter -1 delta 0.00001 --test