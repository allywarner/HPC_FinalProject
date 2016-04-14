objects = matvec.o lanczos.o

lanczos: $(objects)
	cc -o lanczos $(objects) -lm

lanczos.o:
matvec.o:

.PHONY: clean
clean:
	rm lanczos $(objects)
