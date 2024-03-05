gomp-only: main-openmp-only.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c
	gcc -fopenmp -std=c99 main-openmp-only.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c -o gomp-only -lm

gcomplete: main-mpi.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c
	mpicc -fopenmp -std=c99 main-mpi.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c -o gcomplete -lm

iomp-only: main-openmp-only.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c
	icc -qopenmp -std=c99 main-openmp-only.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c -o iomp-only -lm 

icomplete: main-mpi.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c
	mpiicc -qopenmp -std=c99 main-mpi.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c -o icomplete -lm
