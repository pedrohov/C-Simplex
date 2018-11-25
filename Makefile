all:
	gcc -o matrix.o -c matrix.c
	gcc -o simplex.o -c simplex.c
	gcc -o simplex matrix.o simplex.o main.c -lm