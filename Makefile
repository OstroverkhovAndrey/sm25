
all:
	mkdir -p build
	g++ -fopenmp  -c prog.c -o  build/prog.o -I. -g
	g++ -c variant.c -o build/variant.o -I. -g
	g++ -c matrix.c -o build/matrix.o -I. -g
	g++ -fopenmp  build/prog.o build/variant.o build/matrix.o -o build/prog -g
	./build/prog