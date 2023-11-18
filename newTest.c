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

void findFarthestNeighbor(double **distanceMatrix, int numOfCoords, int currentNode, int *farthestIndex, double *farthestDistance) {
    double localFarthestDistance = -DBL_MAX;
    int localFarthestIndex = -1;

    #pragma omp parallel for reduction(max : localFarthestDistance)
    for (int i = 0; i < numOfCoords; i++) {
        if (i != currentNode) {
            double dist = distanceMatrix[currentNode][i];
            if (dist > localFarthestDistance) {
                localFarthestDistance = dist;
            }
        }
    }

    // Find the index associated with the found maximum distance
    for (int i = 0; i < numOfCoords; i++) {
        if (i != currentNode && distanceMatrix[currentNode][i] == localFarthestDistance) {
            localFarthestIndex = i;
            break;
        }
    }

    *farthestDistance = localFarthestDistance;
    *farthestIndex = localFarthestIndex;
}

int *farthestInsertion(double **distanceMatrix, int numOfCoords) {
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

    // Step 2: Find the farthest neighbor to 0
    double farthestDistance;
    int farthestIndex;
    findFarthestNeighbor(distanceMatrix, numOfCoords, 0, &farthestIndex, &farthestDistance);

    tour[1] = farthestIndex;
    tour[2] = 0; // Close the loop
    tourSize = 3;

    while (tourSize < numOfCoords + 1) {
        double globalMaxCost = -DBL_MAX;
        int globalMaxCostIndex = -1, globalInsertPosition = -1;
        #pragma omp parallel
        {
            int threadID = omp_get_thread_num();
            double localMaxCost = -DBL_MAX;
            int localMaxCostIndex = -1, localInsertPosition = -1;

            #pragma omp for nowait
            for (int idx = 0; idx < unvisitedCount; idx++) {
                int i = unvisited[idx];
                for (int j = 0; j < tourSize - 1; j++) {
                    double cost = distanceMatrix[tour[j]][i] + distanceMatrix[i][tour[j + 1]] - distanceMatrix[tour[j]][tour[j + 1]];
                    if (cost > localMaxCost) {
                        localMaxCost = cost;
                        localMaxCostIndex = i;
                        localInsertPosition = j + 1;
                    }
                }
            }

            #pragma omp critical
            {
                if (localMaxCost > globalMaxCost) {
                    globalMaxCost = localMaxCost;
                    globalMaxCostIndex = localMaxCostIndex;
                    globalInsertPosition = localInsertPosition;
                }
            }
        }

        // Update the tour outside the parallel region
        if (globalMaxCostIndex != -1) {
            // Shift elements to the right to make space for the new node
            for (int i = tourSize; i > globalInsertPosition; i--) {
                tour[i] = tour[i - 1];
            }
            // Insert the new node
            tour[globalInsertPosition] = globalMaxCostIndex;
            tourSize++;

            // Remove the visited node from the unvisited array
            for (int i = 0; i < unvisitedCount; i++) {
                if (unvisited[i] == globalMaxCostIndex) {
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
    clock_t start_time = clock();
    // Generate distance matrix
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    // Apply the farthest insertion algorithm
    int *tour = farthestInsertion(distanceMatrix, numOfCoords);

    // Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFile);
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    FILE *file = fopen("timelapsed", "w"); // Open the file once outside the loop
    if (file == NULL) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    fprintf(file, "Elapsed time: %f seconds\n", elapsed_time); // Example of writing to the file

    fclose(file); // Close the file
    // Free memory
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }

    free(distanceMatrix);
    free(tour);

    return 0;
}
