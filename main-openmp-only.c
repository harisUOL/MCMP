#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include<float.h>
#include<math.h>
#include<string.h>
#include<omp.h>

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

	if(argc != 5){
		printf("Program should be called as ./program <coordFile> <outFileC> <outFileF> <outFileN>");
		exit(EXIT_FAILURE);
	}

	//Argument setup for file and output
	char filename[500];
	char outFileC[500];
	char outFileF[500];
	char outFileN[500];

	strcpy(filename, argv[1]);
	strcpy(outFileC, argv[2]);
	strcpy(outFileF, argv[3]);
	strcpy(outFileN, argv[4]);

	//Reading files
	int numOfCoords = readNumOfCoords(filename);
	double **coords = readCoords(filename, numOfCoords);

	//Start program timer
	double tStart = omp_get_wtime();

	/*Program starts*/

	double **dMatrix = createDistanceMatrix(coords, numOfCoords);
	int *currTour;
	int *bestCtour;
	int *bestFtour;
	int *bestNtour;
	double minCcost = __DBL_MAX__;
	double minFcost = __DBL_MAX__;
	double minNcost = __DBL_MAX__;
	double tolerance = 1e-9;

	double tMatrix = omp_get_wtime();
	for(int i =0; i < numOfCoords; i++){
		currTour = solveCheapestInsertion(dMatrix, numOfCoords, i);
		double cost = tourCost(currTour, (numOfCoords+1), dMatrix);
		if(cost+tolerance < minCcost)
		{
			minCcost = cost;
			bestCtour = currTour;
		}
	}

	for(int i =0; i < numOfCoords; i++){
		currTour = solveFarthestInsertion(dMatrix, numOfCoords, i);
		double cost = tourCost(currTour, (numOfCoords+1), dMatrix);
		if(cost+tolerance < minFcost)
		{
			minFcost = cost;
			bestFtour = currTour;
		}
	}

	for(int i =0; i < numOfCoords; i++){
		currTour = solveNearestAddition(dMatrix, numOfCoords, i);
		double cost = tourCost(currTour, (numOfCoords+1), dMatrix);
		if(cost+tolerance < minNcost)
		{
			minNcost = cost;
			bestNtour = currTour;
		}
	}

	double tEnd = omp_get_wtime();
	/*Program ends*/

	printf("\nProgram took %f milliseconds\n", (tEnd - tStart) * 1000);

	if (writeTourToFile(bestCtour, numOfCoords + 1, outFileC) == NULL){
		printf("Error");
	}
	if (writeTourToFile(bestFtour, numOfCoords + 1, outFileF) == NULL){
		printf("Error");
	}
	if (writeTourToFile(bestNtour, numOfCoords + 1, outFileN) == NULL){
		printf("Error");
	}

	//Free memory
	for(int i = 0; i < numOfCoords; i++){
		free(dMatrix[i]);
	}

	free(dMatrix);
}