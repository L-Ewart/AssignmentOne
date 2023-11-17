#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>

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
        for (int j = 0; j < numOfCoords; j++) {
            matrix[i][j] = distance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
        }
    }
    return matrix;
}

int *furthestInsertion(double **distanceMatrix, int numOfCoords) {
    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int));
    bool *visited = (bool *)calloc(numOfCoords, sizeof(bool));
    int tourSize = 1;

    // Step 1: Start with vertex 0
    tour[0] = 0;
    visited[0] = true;

    // Step 2: Find the furthest neighbor to 0
    int furthest = 0, furthestIndex = -1;
    for (int i = 1; i < numOfCoords; i++) {
        if (distanceMatrix[0][i] > furthest) {
            furthest = distanceMatrix[0][i];
            furthestIndex = i;
        }
    }
    tour[1] = furthestIndex;
    tour[2] = 0; // Close the loop
    visited[furthestIndex] = true;
    tourSize = 3;

    // Step 3 and 4: Repeat until all vertices are visited
    while (tourSize < numOfCoords + 1) {
        int maxCost = 0;
        int maxCostIndex = -1, insertPosition = -1;

        for (int i = 0; i < numOfCoords; i++) {
            if (!visited[i]) {
                for (int j = 0; j < tourSize - 1; j++) {
                    int cost = distanceMatrix[tour[j]][i] + distanceMatrix[i][tour[j + 1]] - distanceMatrix[tour[j]][tour[j + 1]];
                    if (cost > maxCost) {
                        maxCost = cost;
                        maxCostIndex = i;
                        insertPosition = j + 1;
                    }
                }
            }
        }

        // Insert the node
        for (int i = tourSize; i > insertPosition; i--) {
            tour[i] = tour[i - 1];
        }
        tour[insertPosition] = maxCostIndex;
        visited[maxCostIndex] = true;
        tourSize++;
    }

    free(visited);
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
    int *tour = furthestInsertion(distanceMatrix, numOfCoords);

    // Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFile);

    // Free memory
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);
    }
    free(coords);
    free(distanceMatrix);
    free(tour);

    return 0;
}