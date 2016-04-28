objects = matvec.o eig.o scan.o partition.o genDot.o quicksort.o

partition: $(objects)
	mpicc -o partition $(objects) -lm -lblas -llapack -fopenmp

partition.o: partition.c
	mpicc -c partition.c -fopenmp
# lanczos.o: lanczos.c
# 	mpicc -c lanczos.c -fopenmp
matvec.o: matvec.c
	cc -c matvec.c -fopenmp -fopenmp
eig.o: eig.c
	cc -c eig.c -llapack -lblas
scan.o: scan.c
	cc -c scan.c -fopenmp
genDot.o: genDot.c
	cc -c genDot.c -fopenmp
quicksort.o: Matrices/quicksort.c
	cc -c Matrices/quicksort.c -fopenmp

.PHONY: clean
clean:
	rm partition $(objects)
