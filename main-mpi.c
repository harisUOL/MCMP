#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include<float.h>
#include<math.h>
#include<string.h>
#include<omp.h>
#include<mpi.h>

// Explicit declaration of functions from "coordReader.c" file
int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);
int *solveCheapestInsertion(double **d, int numOfCoords, int start);
int *solveFarthestInsertion(double **d, int numOfCoords, int start);
int *solveNearestAddition(double **d, int numOfCoords, int start);

double tourCost(int *tour, int tourLength, double **distance_matrix){
	double cost = 0;

	for(int i = 1; i < tourLength; i++){
		cost += distance_matrix[tour[i - 1]][tour[i]];
	}
	return cost;
}


double **createDistanceMatrix(double **coords, int numOfCoords){
	int i, j;
	
	double **dMatrix = (double **)malloc(numOfCoords * sizeof(double));

	for(i = 0; i < numOfCoords; i++){
		dMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
	}

	#pragma omp parallel for collapse(2)
	for(i = 0; i < numOfCoords; i++){
		for(j = 0; j < numOfCoords; j++){
			double diffX = coords[i][0] - coords[j][0];
			double diffY = coords[i][1] - coords[j][1];
			dMatrix[i][j] = sqrt((diffX * diffX) + (diffY * diffY));
		}
	}

	return dMatrix;
}

int main(int argc, char *argv[]){

	int myRank, nProcs;

	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	int bCastRoot = 0;
	int numOfCoords, myStart, myEnd, *bestCtours, *bestFtours, *bestNtours;
	int allStarts[nProcs], allEnds[nProcs];
	char filename[500];
	char outFileC[500];
	char outFileF[500];
	char outFileN[500];
	double tStart, tEnd;
	double flatDistMatrix[512*512];
	double **dMatrix;

	if(myRank == bCastRoot)
	{
		if(argc != 5){
			printf("Program should be called as ./program <coordFile> <outFileC> <outFileF> <outFileN>");
			exit(EXIT_FAILURE);
		}

		//Argument setup for file and output

		strcpy(filename, argv[1]);
		strcpy(outFileC, argv[2]);
		strcpy(outFileF, argv[3]);
		strcpy(outFileN, argv[4]);

		//Reading files
		numOfCoords = readNumOfCoords(filename);
		double **coords = readCoords(filename, numOfCoords);

		//Start program timer
		tStart = MPI_Wtime();

		int block = numOfCoords/nProcs;
		int rem = numOfCoords%nProcs;
		for (int i = 0; i < nProcs; ++i)
		{
			allStarts[i] = i*block;
			allEnds[i] = (i+1)*block - 1;
			// Let last proc do the remainder work
			if(i==(nProcs-1))
				allEnds[i] += rem;
		}
		dMatrix = createDistanceMatrix(coords, numOfCoords);
		for (int i = 0; i < numOfCoords; ++i)
		{
			for (int j = 0; j < i; ++j)
				flatDistMatrix[j+(numOfCoords*i)] = flatDistMatrix[i+(numOfCoords*j)] = dMatrix[i][j];
			flatDistMatrix[i+(numOfCoords*i)] = 0;
		}
		bestCtours = (int *) malloc((numOfCoords+1)*nProcs*sizeof(int));
		bestFtours = (int *) malloc((numOfCoords+1)*nProcs*sizeof(int));
		bestNtours = (int *) malloc((numOfCoords+1)*nProcs*sizeof(int));
	}
	/*Program starts*/

	MPI_Bcast(&numOfCoords, 1, MPI_INT, bCastRoot, MPI_COMM_WORLD);
	MPI_Bcast(&allStarts, nProcs, MPI_INT, bCastRoot, MPI_COMM_WORLD);
	MPI_Bcast(&allEnds, nProcs, MPI_INT, bCastRoot, MPI_COMM_WORLD);
	MPI_Bcast(&flatDistMatrix, numOfCoords*numOfCoords, MPI_DOUBLE, bCastRoot, MPI_COMM_WORLD);

	// Unpack flat matrix
	dMatrix = (double **)malloc(numOfCoords*sizeof(double));
	for (int i = 0; i < numOfCoords; ++i)
		dMatrix[i] = (double *)malloc(numOfCoords*sizeof(double));

	for (int i = 0; i < numOfCoords; ++i)
	{
		for (int j = 0; j < numOfCoords; ++j)
			dMatrix[i][j] = flatDistMatrix[i + (j*numOfCoords)];
	}

	myStart = allStarts[myRank];
	myEnd = allEnds[myRank];

	int *currTour, *bestCtour, *bestFtour, *bestNtour;
	double minCcost = __DBL_MAX__;
	double minFcost = __DBL_MAX__;
	double minNcost = __DBL_MAX__;
	double tolerance = 1e-9;

	for(int i = myStart; i <= myEnd; i++){
		currTour = solveCheapestInsertion(dMatrix, numOfCoords, i);
		double cost = tourCost(currTour, (numOfCoords+1), dMatrix);
		if(cost+tolerance < minCcost)
		{
			minCcost = cost;
			bestCtour = currTour;
		}
	}

	for(int i = myStart; i <= myEnd; i++){
		currTour = solveFarthestInsertion(dMatrix, numOfCoords, i);
		double cost = tourCost(currTour, (numOfCoords+1), dMatrix);
		if(cost+tolerance < minFcost)
		{
			minFcost = cost;
			bestFtour = currTour;
		}
	}

	for(int i = myStart; i < myEnd; i++){
		currTour = solveNearestAddition(dMatrix, numOfCoords, i);
		double cost = tourCost(currTour, (numOfCoords+1), dMatrix);
		if(cost+tolerance < minNcost)
		{
			minNcost = cost;
			bestNtour = currTour;
		}
	}

	// Gather answers
	MPI_Gather(bestCtour, (numOfCoords+1), MPI_INT, bestCtours, (numOfCoords+1), MPI_INT, bCastRoot, MPI_COMM_WORLD);
	MPI_Gather(bestFtour, (numOfCoords+1), MPI_INT, bestFtours, (numOfCoords+1), MPI_INT, bCastRoot, MPI_COMM_WORLD);
	MPI_Gather(bestNtour, (numOfCoords+1), MPI_INT, bestNtours, (numOfCoords+1), MPI_INT, bCastRoot, MPI_COMM_WORLD);

	if(myRank == bCastRoot)
	{
		int bestCIn, bestFIn, bestNIn;
		int *currTourC = malloc((numOfCoords+1)*sizeof(int));
		int *currTourF = malloc((numOfCoords+1)*sizeof(int));
		int *currTourN = malloc((numOfCoords+1)*sizeof(int));
		minCcost = minFcost = minNcost = __DBL_MAX__;
		double costC,costF,costN;
		for (int i = 0; i < nProcs; ++i)
		{
			for (int j = 0; j <= numOfCoords; ++j)
			{
				currTourC[j] = bestCtours[i*(numOfCoords+1) + j];
				currTourF[j] = bestFtours[i*(numOfCoords+1) + j];
				currTourN[j] = bestNtours[i*(numOfCoords+1) + j];
			}
			costC = tourCost(currTourC, (numOfCoords+1), dMatrix);
			costF = tourCost(currTourF, (numOfCoords+1), dMatrix);
			costN = tourCost(currTourN, (numOfCoords+1), dMatrix);
			if( (costC + tolerance) < minCcost )
			{
				minCcost = costC;
				bestCIn = i;
			}
			if( (costF + tolerance) < minFcost )
			{
				minFcost = costF;
				bestFIn = i;
			}
			if( (costN + tolerance) < minNcost )
			{
				minNcost = costN;
				bestNIn = i;
			}
		}
		for (int j = 0; j <= numOfCoords; ++j)
		{
			currTourC[j] = bestCtours[bestCIn*(numOfCoords+1) + j];
			currTourF[j] = bestFtours[bestFIn*(numOfCoords+1) + j];
			currTourN[j] = bestNtours[bestNIn*(numOfCoords+1) + j];
		}

		if (writeTourToFile(currTourC, numOfCoords + 1, outFileC) == NULL){
			printf("Error");
		}
		if (writeTourToFile(currTourF, numOfCoords + 1, outFileF) == NULL){
			printf("Error");
		}
		if (writeTourToFile(currTourN, numOfCoords + 1, outFileN) == NULL){
			printf("Error");
		}
		double tEnd = MPI_Wtime();
		/*Program ends*/

		printf("\nProgram took %f milliseconds\n", (tEnd - tStart) * 1000);

		free(bestCtours);
		free(bestFtours);
		free(bestNtours);
		free(currTourC);
		free(currTourF);
		free(currTourN);
	}


	//Free memory
	for(int i = 0; i < numOfCoords; i++){
		free(dMatrix[i]);
	}

	free(dMatrix);
}