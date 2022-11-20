FLAGS = -g

all: xsimpl test_Simplex 

test_Simplex: nrutil.o simp1.o simp2.o simp3.o simplx.o Simplex.o test_Simplex.cpp
	g++ ${FLAGS} nrutil.o simp1.o simp2.o simp3.o simplx.o Simplex.o test_Simplex.cpp -o test_Simplex

Simplex.o: Simplex.cpp Simplex.h nrutil.o simp1.o simp2.o simp3.o simplx.o
	g++ -c Simplex.cpp

xsimpl: nrutil.o simp1.o simp2.o simp3.o simplx.o xsimplx.c
	gcc  -lm nrutil.o  simp1.o simp2.o simp3.o simplx.o xsimplx.c  -o xsimpl

simplx.o: simplx.c 
	gcc -lm -c simplx.c

nrutil.o: nrutil.c 
	gcc -lm -c nrutil.c

simp1.o: simp1.c
	gcc -lm -c simp1.c

simp2.o: simp2.c
	gcc -lm -c simp2.c

simp3.o: simp3.c
	gcc -lm -c simp3.c



clean:
	rm -f *.o xsimpl test_extern test_Simplex a.out a

run_all: all
	echo "\nSimplex de C:"
	./xsimpl
	echo "\nSimplex de C++:"
	./test_Simplex