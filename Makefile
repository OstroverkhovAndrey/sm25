
all:
	mkdir -p build
	g++ -fopenmp  -c prog.c -o  build/prog.o -I. -g
	g++ -c variant.c -o build/variant.o -I. -g
	g++ -c matrix.c -o build/matrix.o -I. -g
	g++ -c args.c -o build/args.o -I. -g
	g++ -fopenmp  build/prog.o build/variant.o build/matrix.o build/args.o -o build/prog -g
	./build/prog N 40 M 60 count_iter -1 delta 0.00001 --test