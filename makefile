objects = matvec.o lanczos.o eig.o

lanczos: $(objects)
	cc -o lanczos $(objects) -lm -lblas -llapack

lanczos.o: lanczos.c
	cc -c lanczos.c -fopenmp
matvec.o: matvec.c
	cc -c matvec.c -fopenmp
eig.o: eig.c
	cc -c eig.c -llapack -lblas

.PHONY: clean
clean:
	rm lanczos $(objects)
