# Parallel and Distributed programming in HPC
<i><b>
This repository consist of three parallelised methods for solving the travelling salesman problem 
which are executed upon multiple number of co-ordinates(being 9, 16, and 512 in numbers) and then 
using distrubuted programming those programs are executed on multiple nodes and with multiple instances instantiating
every co-ordinate within the co-ordinates files as the first and last points of the TSP tours. Submitting these files through HPC makes it 
quicker and easier to run those many instances at once and through this we find the best tour out of all the
returned tours.
</i></b>
<hr>
<br>

### # Files included within are as follows:
- coordReader.c: includes fucntions that read the coordinate input from a file, and write the final tour to a file.
- 9_coords.coord, 16_coords.coord and 512_coords.coord: co-ordinate files.
- 9_cout.dat, 9_fout.dat, 9_nout.dat, 16_cout.dat, ... : expected output files from the three methods(cInsertion, fInsertion, and nAddition).
- Continuous Assessment files 1 and 2: files including outline of the parallelism and distribution of methods.
- ompcInsertion.c, ompfInsertion.c, ompnAddition.c: parallelised versions of methods implemented in c
- main-openmp-only.c: sequential implementation of all three methods to give out the best tour starting from one co-cordinate only.
- main-mpi.c: distributed implementation of the parallelised program running those methods on multiple instances through multiple nodes.
- Makefile: compiles all versions of c and makes them into executables to be ran in terminals.
- MPI_batch.sh: bash script to submit jobs on Barkla server(in this case, University of Liverpool's Server) through slurm.
- Output folder: consists of output files given through the bash script which shows the execution of programs on different coordinates and there execution time and the number of threads used.

### # Speedup plot
<br>

![Screenshot from 2024-01-20 21-46-30](https://github.com/harisUOL/MCMP/assets/146499851/887f77d1-6fe5-4147-a073-57314c8ed241)

### # Parallel efficiency plot
<br>

![Screenshot from 2024-01-20 22-16-03](https://github.com/harisUOL/MCMP/assets/146499851/649c01c4-1360-4d6f-a9ff-48c3c0119003)
