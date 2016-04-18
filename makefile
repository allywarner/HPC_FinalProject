objects = matvec.o lanczos.o eig.o scan.o

lanczos: $(objects)
	cc -o lanczos $(objects) -lm -lblas -llapack -fopenmp

lanczos.o: lanczos.c
	cc -c lanczos.c -fopenmp
matvec.o: matvec.c
	cc -c matvec.c -fopenmp -fopenmp
eig.o: eig.c
	cc -c eig.c -llapack -lblas
scan.o: scan.c
	cc -c scan.c -fopenmp

.PHONY: clean
clean:
	rm lanczos $(objects)
