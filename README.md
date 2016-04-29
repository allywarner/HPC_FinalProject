## Generate P Matrix
To run the program, you first need to generate a p matrix with the preprocessor in the Matrices folder.  All you need to do is download the mtx file from the website, then remove all the stuff at the top above the row that lists the dimensions and number of lines in the file.  Then just run ./preprocessor <path-to-matrix>.  This will generate a file with the same name but a p_ in front.

## partition the graph

Then all you have to do is log into chpc and run
sbatch slurm.sh

The very last line of slurm.sh is

mpirun ./partition Matrices/p_finance256.dat 100 1

The 100 is just the number of iterations you want in the lanczos algorithm.  The 1 is just whether or not you want to do full reorthogonalization.

