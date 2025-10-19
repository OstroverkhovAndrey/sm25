
all:
	g++ -fopenmp  -c proc.c -o  proc.o -I. -g
	g++ -c variant.c -o variant.o -I. -g
	g++ -c matrix.c -o matrix.o -I. -g
	g++ -fopenmp  proc.o variant.o matrix.o -o proc -g
	./proc