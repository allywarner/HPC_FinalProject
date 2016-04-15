objects = matvec.o lanczos.o eig.o

lanczos: $(objects)
	cc -o lanczos $(objects) -lm -lblas -llapack

lanczos.o:
matvec.o:
eig.o: eig.c
	cc -c eig.c -llapack -lblas

.PHONY: clean
clean:
	rm lanczos $(objects)
