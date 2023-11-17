#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <omp.h>

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

int min_distance = INT_MAX;

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



//#define N 4 // Number of cities

//int graph[N][N] = {
    //{0, 10, 15, 20},
    //{10, 0, 35, 25},
    //{15, 35, 0, 30},
   // {20, 25, 30, 0}
//};


int order[numOfCoords];

void TSP(int start_city, int numOfCoord) {
    int visited[numOfCoord] = {0};
    int current_city = start_city;
    visited[start_city] = 1;

    int count = 1;
    int distance = 0;
    int path[numOfCoord];
    path[0] = start_city;

    while (count < numOfCoord) {
        int nearest = INT_MAX;
        int next_city = -1;

        #pragma omp parallel for
        for (int i = 0; i < numOfCoord; i++) {
            if (graph[current_city][i] < nearest && !visited[i]) {
                #pragma omp critical
                {
                    if (graph[current_city][i] < nearest) {
                        nearest = graph[current_city][i];
                        next_city = i;
                    }
                }
            }
        }

        distance += nearest;
        visited[next_city] = 1;
        current_city = next_city;
        path[count] = next_city;
        count++;
    }

    distance += graph[current_city][start_city]; // Complete the cycle

    #pragma omp critical
    {
        if (distance < min_distance) {
            min_distance = distance;
            for (int i = 0; i < N; i++) {
                order[i] = path[i];
            }
        }
    }
}

int main(int argc, char *argv[]) {

    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    char *inputFile = argv[1];
    char *outputFile = argv[2];

    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);

    // Generate distance matrix
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);


    #pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        TSP(i, numOfCoords);
    }

    printf("Minimum distance: %d\n", min_distance);
    printf("Order of cities visited: ");
    for (int i = 0; i < numOfCoords; i++) {
        printf("%d ", order[i]);
    }
    printf("\n");

    return 0;
}