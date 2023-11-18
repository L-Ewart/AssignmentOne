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

// Function to find the nearest neighbor
void findNearestNeighbor(double **distanceMatrix, int numOfCoords, int currentNode, int *nearestIndex, double *nearestDistance) {
    double localNearestDistance = DBL_MAX;
    int localNearestIndex = -1;

    #pragma omp parallel for reduction(min : localNearestDistance)
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

    tour[1] = nearestIndex;
    tour[2] = 0; // Close the loop
    tourSize = 3;

    while (tourSize < numOfCoords + 1) {
        double globalMinCost = DBL_MAX;
        int globalMinCostIndex = -1, globalInsertPosition = -1;

        #pragma omp parallel
        {
            int threadID = omp_get_thread_num();
            double localMinCost = DBL_MAX;
            int localMinCostIndex = -1, localInsertPosition = -1;

            #pragma omp for nowait
            for (int idx = 0; idx < unvisitedCount; idx++) {
                int i = unvisited[idx];
                for (int j = 0; j < tourSize - 1; j++) {
                    double cost = distanceMatrix[tour[j]][i] + distanceMatrix[i][tour[j + 1]] - distanceMatrix[tour[j]][tour[j + 1]];
                    if (cost < localMinCost) {
                        localMinCost = cost;
                        localMinCostIndex = i;
                        localInsertPosition = j + 1;
                    }
                }
            }

            #pragma omp critical
            {
                if (localMinCost < globalMinCost) {
                    globalMinCost = localMinCost;
                    globalMinCostIndex = localMinCostIndex;
                    globalInsertPosition = localInsertPosition;
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
                    unvisitedCount--;                               // Reduce the count of unvisited nodes
                    break;
                }
            }
        }
    }

    free(unvisited);
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

    // Generate distance matrix
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    // Apply the cheapest insertion algorithm
    int *tour = cheapestInsertion(distanceMatrix, numOfCoords);

    // Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFile);

    // Free memory
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
    free(tour);

    return 0;
}
