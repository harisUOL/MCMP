#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>
#include<string.h>
#include<omp.h>

int *solveCheapestInsertion(double **dMatrix, int numOfCoords, int start){

	//Memory allocation for the tour and visited arrays. Tour is numOfCoords + 1 for returning to origin
	int *tour = (int *)malloc((1 + numOfCoords) * sizeof(int));
	//Visited uses calloc, array is instantiated with "0" as all elements. Good for boolean arrays.
	bool *visited = (bool *)calloc(numOfCoords, sizeof(bool));

	//Initialising tour to empty
	for(int i = 0; i < numOfCoords; i++){
		tour[i] = -1;
	}

	//Tour always starts and ends at 0. 0 is visited
	tour[0] = start;
	tour[1] = start;
	visited[start] = true;
	
	//Hard coding because I'm lazy
	int numVisited = 1;
	int tourLength = 2;

	//Where OMP starts... Get the env variable for the max num of threads.
	int numThreads = omp_get_max_threads();

	double *threadMinCosts = NULL;
	int *threadNextNode = NULL;
	int *threadInsertPos = NULL;
		
	threadMinCosts = (double*)malloc(numThreads * 8 * sizeof(double));
	threadNextNode = (int*)malloc(numThreads * 8 * sizeof(int));
	threadInsertPos = (int*)malloc(numThreads * 8 * sizeof(int));

	//Start a parallel section
	#pragma omp parallel 
	{
		//Each thread now has started, and it stores its thread number in threadID
		int threadID = omp_get_thread_num();
		while(numVisited < numOfCoords){

			//Thread only accesses its memory location in the shared array. No race conditions.
			threadMinCosts[threadID * 8] = __DBL_MAX__;
			threadNextNode[threadID * 8] = -1;
			threadInsertPos[threadID * 8] = -1;

			//Begin a workshare construct. Threads divide i and j and work on their respective iterations.
			#pragma omp for collapse(2)
			for(int i = 0; i < tourLength - 1; i++){	
				for(int j = 0; j < numOfCoords; j++){

					//Each thread performs their cheapest insertion. Works on each position in the tour.
					if(!visited[j]){
						double cost = dMatrix[tour[i]][j] + dMatrix[tour[i+1]][j] - dMatrix[tour[i]][tour[i + 1]];
						if(cost < threadMinCosts[threadID * 8]){

							threadMinCosts[threadID * 8] = cost;
							threadNextNode[threadID * 8] = j;
							threadInsertPos[threadID * 8] = i + 1;
						}
					}
				}
			}

			//Only one thread works on this part. This part must be serial. OMP single instead of master. Therefore implicit barrier
			#pragma omp single
			{
				int bestNextNode = -1;
				int bestInsertPos = -1;
				double minCost = __DBL_MAX__;

				//A single thread loops through each threads memory locations. Finds the minCost
				for(int i = 0; i < numThreads; i++){
					if(threadMinCosts[i * 8] < minCost){
						minCost = threadMinCosts[i * 8];
						bestNextNode = threadNextNode[i * 8];
						bestInsertPos = threadInsertPos[i * 8];
					}
				}	

				//One thread places the bestNextNode in the bestInsertPos
				for(int i = numOfCoords; i > bestInsertPos; i--){
					tour[i] = tour[i - 1];
				}

				tour[bestInsertPos] = bestNextNode;
				visited[bestNextNode] = true;		
				
				tourLength++;
				numVisited++;
			}
		}
	}

	//Free all memory when done
	
	free(visited);
	free(threadMinCosts);
	free(threadNextNode);
	free(threadInsertPos);

	return tour;
}