#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>

// Function declarations
int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);

int find_index(int arr[], int size, int element) {
    for (int i = 0; i < size; i++) {
        if (arr[i] == element) {
            return i;
        }
    }
    return -1; // element not found
}
void calculateDistanceMatrix(double **coords, int numOfCoords, double **distM);

int main(int argc, char *argv[]) {
    // Check for the correct number of command-line arguments
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file> <output_file>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    char *filename = argv[1];
    char *outputFile = argv[2];

    int numOfCoords = readNumOfCoords(filename);
    printf("Number of Coordinates: %d\n", numOfCoords);

    // Allocate memory for coordinates and distance matrix
    double **coords = readCoords(filename, numOfCoords);
    double **distM = (double **)malloc(numOfCoords * sizeof(double *));
    for (int j = 0; j < numOfCoords; j++) {
        distM[j] = (double *)malloc(numOfCoords * sizeof(double));
    }

    // Calculate the distance matrix
    calculateDistanceMatrix(coords, numOfCoords, distM);

    // Array to keep track of the available options
    int *toVisit = (int *)malloc(numOfCoords * sizeof(int));
    for (int counter = 0; counter < numOfCoords; counter++) {
        toVisit[counter] = counter;
    }

    // Array to store the steps taken during the journey
    int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int));
    tour[0] = 0;

    // Main loop for visiting each location
    for (int visitNumber = 0; visitNumber < numOfCoords - 1; visitNumber++) {
        int minPosition;
        double minCost = INFINITY;

        // Loop to find the minimum cost and position to visit
        #pragma omp parallel for shared(minPosition, minCost)
        for (int nextCheck = 0; nextCheck < numOfCoords; nextCheck++) {
            int nextPosition = toVisit[nextCheck];

            if (nextPosition != 0) {
                for (int positionBefore = 0; positionBefore < visitNumber + 1; positionBefore++) {
                    double minimalCost = distM[tour[positionBefore]][nextPosition] +
                                         distM[nextPosition][tour[positionBefore + 1]] -
                                         distM[tour[positionBefore]][tour[positionBefore + 1]];

                    #pragma omp critical
                    {
                        if (minimalCost < minCost) {
                            minCost = minimalCost;
                            minPosition = nextPosition;
                        }
                    }
                }
            }
        }

        // Update the tour based on the minimum cost and position
        int indexA = find_index(tour, visitNumber + 2, tour[visitNumber]);
        int indexB = (indexA == visitNumber + 1) ? indexA : indexA + 1;

        for (int i = visitNumber; i >= indexB; i--) {
            int temp = tour[i];
            tour[i + 1] = temp;
        }

        tour[indexA + 1] = minPosition;
        toVisit[minPosition] = 0;

        // Print the visiting order
        printf("Visiting Order: ");
        for (int i = 0; i < numOfCoords + 1; i++) {
            printf("%d ", tour[i]);
        }
        printf("\n");
    }

    // Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFile);

    // Free allocated memory
    free(toVisit);
    free(tour);

    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
        free(distM[i]);
    }
    free(coords);
    free(distM);

    return 0;
}

// Function to calculate the distance matrix
void calculateDistanceMatrix(double **coords, int numOfCoords, double **distM) {
    #pragma omp for
    for (int xx = 0; xx < numOfCoords; xx++) {
        for (int yy = 0; yy < numOfCoords; yy++) {
            distM[xx][yy] = sqrt(pow(coords[yy][0] - coords[xx][0], 2) + pow(coords[yy][1] - coords[xx][1], 2));
        }
    }

    // Uncomment to print the distance matrix
    printf("Distance Matrix\n");
    for (int x = 0; x < numOfCoords; x++) {
        for (int y = 0; y < numOfCoords; y++) {
            printf("%f    ", distM[x][y]);
        }
        printf("\n");
    }
}