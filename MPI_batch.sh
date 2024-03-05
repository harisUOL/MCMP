#!/bin/bash -l

# Specific course queue and max wallclock time
#SBATCH -p course -t 15

# Defaults on Barkla (but set to be safe)
## Specify the current working directory as the location for executables/files
#SBATCH -D ./
## Export the current environment to the compute node
#SBATCH --export=ALL

# load modules intel compiler
module load compilers/intel/2019u5 
module load mpi/intel-mpi/2019u5/bin

echo "Node list                    : $SLURM_JOB_NODELIST"
echo "Number of nodes allocated    : $SLURM_JOB_NUM_NODES or $SLURM_NNODES"
echo "Number of threads or processes          : $SLURM_NTASKS"
echo "Number of processes per node : $SLURM_TASKS_PER_NODE"
echo "Requested tasks per node     : $SLURM_NTASKS_PER_NODE"
echo "Requested CPUs per task      : $SLURM_CPUS_PER_TASK"
echo "Scheduling priority          : $SLURM_PRIO_PROCESS"

SRC=$1
A1=$2
A2=$3
A3=$4
COMM=$5

#sets the exe name as the sourcecode, and %% removes the ".c"
EXE=gcomplete
PROCS=${SLURM_NTASKS:-1}

#sets test
TEST1=$6_coords.coord
OUT1=$6_output_c_mpi.txt
OUT2=$6_output_f_mpi.txt
OUT3=$6_output_n_mpi.txt

#deletes the existing executable (if it exists)
rm -f ${EXE}

echo compiling $SRC to $EXE 

#compilation using intel compiler of sourcecode to exectuable.
mpicc $SRC $A1 $A2 $A3 $COMM -o $EXE -std=c99 -lm -fopenmp
echo
echo ------------------------------------

#this tests if the file $EXE is present and executable
if test -x $EXE; then 
      # set number of threads
      
      # if '-c' not used then default to 1. SLURM_CPUS_PER_TASK is given by -c
      export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} 
      echo using ${PROCS} processes
      echo using ${OMP_NUM_THREADS} OpenMP threads
      echo 
      echo Multiple execution..
      echo
      echo

      # run multiple times. Because we have exported how many threads we're using, we just execute the file.
      for i in {1..1}; do mpirun ./${EXE} $TEST1 $OUT1 $OUT2 $OUT3; done
else
     echo $SRC did not built to $EXE
fi
