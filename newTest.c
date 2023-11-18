#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <omp.h>
#include <float.h>
#include <time.h>

#ifndef coordReader_H_
#define coordReader_H_

int readNumOfCoords();
double **readCoords();
void *writeTourToFile();
#endif
// Provided functions
int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);


// Custom reduction function
/* MinCostData minCostFunc(MinCostData a, MinCostData b) {
    return (a.cost < b.cost) ? a : b;
} */


// Function to calculate the Euclidean distance between two points
double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

// Function to generate a distance matrix from coordinates



double **generateDistanceMatrix(double **coords, int numOfCoords) {
    double **matrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }
    #pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i][i] = 0; // Distance from a point to itself is 0
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix[i][j] = distance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix[j][i] = matrix[i][j]; // Use symmetry, avoid redundant calculation
        }
    }
    return matrix;
}

     void findNearestNeighbor(double **distanceMatrix, int numOfCoords, int currentNode, int *nearestIndex, double *nearestDistance) {
    double localNearestDistance = DBL_MAX;
    int localNearestIndex = -1;

    // Parallel for loop to find the minimum distance
    #pragma omp parallel for reduction(min:localNearestDistance)
    for (int i = 0; i < numOfCoords; i++) {
        if (i != currentNode) {
            double dist = distanceMatrix[currentNode][i];
            if (dist < localNearestDistance) {
                localNearestDistance = dist;
            }
        }
    }

    // Find the index associated with the found minimum distance
    for (int i = 0; i < numOfCoords; i++) {
        if (i != currentNode && distanceMatrix[currentNode][i] == localNearestDistance) {
            localNearestIndex = i;
            break;
        }
    }

    *nearestDistance = localNearestDistance;
    *nearestIndex = localNearestIndex;
} 

int *cheapestInsertion(double **distanceMatrix, int numOfCoords) {
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    int *unvisited = malloc(numOfCoords * sizeof(int));
    int unvisitedCount = numOfCoords - 1;
    
    // Initialize unvisited nodes
    for (int i = 1; i < numOfCoords; i++) {
        unvisited[i - 1] = i;
    }
    // Step 1: Start with vertex 0
    tour[0] = 0;
    int tourSize = 1;
    // Step 2: Find the nearest neighbor to 0
    double nearestDistance;
    int nearestIndex;
    findNearestNeighbor(distanceMatrix, numOfCoords, 0, &nearestIndex, &nearestDistance);  
        //int nearest = INT_MAX, nearestIndex = -1;
    tour[1] = nearestIndex;
    tour[2] = 0; // Close the loop
    tourSize = 3;
    

    int numThreads = omp_get_max_threads();
    double *localMinCosts = malloc(numThreads * sizeof(double));
    int *localMinCostIndices = malloc(numThreads * sizeof(int));
    int *localInsertPositions = malloc(numThreads * sizeof(int));
    while (tourSize < numOfCoords + 1) {
        double globalMinCost = DBL_MAX; 
        int globalMinCostIndex = -1, globalInsertPosition = -1;
        #pragma omp parallel
        {

            int threadID = omp_get_thread_num();
            localMinCosts[threadID] = DBL_MAX;
            localMinCostIndices[threadID] = -1;
            localInsertPositions[threadID] = -1;

            #pragma omp for nowait
            for (int idx = 0; idx < unvisitedCount; idx++) {
                        int i = unvisited[idx];
                    for (int j = 0; j < tourSize - 1; j++) {
                        int cost = distanceMatrix[tour[j]][i] + distanceMatrix[i][tour[j + 1]] - distanceMatrix[tour[j]][tour[j + 1]];
                        if (cost < localMinCosts[threadID]) {
                            localMinCosts[threadID] = cost;
                            localMinCostIndices[threadID] = i;
                            localInsertPositions[threadID] = j + 1;
                        }
                    }
                
            }
            #pragma omp critical
            {
                if (localMinCosts[threadID] < globalMinCost) {
                globalMinCost = localMinCosts[threadID];
                globalMinCostIndex = localMinCostIndices[threadID];
                globalInsertPosition = localInsertPositions[threadID];
                }
            }
        }
        
           
        
        
    

                // Update the tour outside the parallel region
            if (globalMinCostIndex != -1) {
            // Shift elements to the right to make space for the new node
            for (int i = tourSize; i > globalInsertPosition; i--) {
                tour[i] = tour[i - 1];
            }
            // Insert the new node
            tour[globalInsertPosition] = globalMinCostIndex;
            tourSize++;

            // Remove the visited node from the unvisited array
            for (int i = 0; i < unvisitedCount; i++) {
                if (unvisited[i] == globalMinCostIndex) {
                    unvisited[i] = unvisited[unvisitedCount - 1]; // Swap with the last element
                    unvisitedCount--; // Reduce the count of unvisited nodes
                    break;
                }
            }
        }
            
        
    }
    
    free(unvisited);
    free(localMinCosts);
    free(localMinCostIndices);
    free(localInsertPositions);
    return tour;
}

 



int main(int argc, char *argv[]) {
    // Ensure correct usage
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    char *inputFile = argv[1];
    char *outputFile = argv[2];

    // Read coordinates
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    clock_t start_time = clock();
    // Generate distance matrix
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    // Apply the cheapest insertion algorithm
    int *tour = cheapestInsertion(distanceMatrix, numOfCoords);

    // Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFile);
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

   FILE *file = fopen("timelapsed", "w"); // Open the file once outside the loop
    if (file == NULL) {
        perror("Error opening file");
    return EXIT_FAILURE;
    }

// Your processing logic here...

fprintf(file, "Elapsed time: %f seconds\n", elapsed_time); // Example of writing to the file

fclose(file); // Close the file
    // Free memory
    for (int i = 0; i < numOfCoords; i++) {
        //free(coords[i]);
        free(distanceMatrix[i]);

    
    }
    
    
    //free(coords);
    free(distanceMatrix);
    free(tour);

    return 0;
}
