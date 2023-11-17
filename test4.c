#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef coordReader_H_
#define coordReader_H_

int readNumOfCoords();
double **readCoords ();
void *writeTourToFile();
#endif

int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
// Function to calculate the Euclidean distance between two points
double calculate_distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

// Function to convert a 2D array of coordinates into a distance matrix
void convert_to_distance_matrix(double coords[][2], int num_coords, double **distance_matrix) {
    for (int i = 0; i < num_coords; i++) {
        for (int j = 0; j < num_coords; j++) {
            // Calculate the distance between points i and j
            distance_matrix[i][j] = calculate_distance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
        }
    }
}

// Function to print the distance matrix
void print_distance_matrix(int num_coords, double **distance_matrix) {
    printf("Distance Matrix:\n");
    for (int i = 0; i < num_coords; i++) {
        for (int j = 0; j < num_coords; j++) {
            printf("%lf\t", distance_matrix[i][j]);
        }
        printf("\n");
    }
}

int main() {
    // Example: 2D array of coordinates (x, y)
    char input_filename[] = "9_coords.coord";
    int num_coords = readNumOfCoords(input_filename);
    printf("Number of coordinates: %d\n", num_coords);

    double coords[][2] = readCoords(input_filename, num_coords);

    // Number of coordinates
    //int num_coords = sizeof(coords) / sizeof(coords[0]);

    // Allocate memory for the distance matrix
    double **distance_matrix = malloc(num_coords * sizeof(double *));
    for (int i = 0; i < num_coords; i++) {
        distance_matrix[i] = malloc(num_coords * sizeof(double));
    }

    // Convert coordinates to distance matrix
    convert_to_distance_matrix(coords, num_coords, distance_matrix);

    // Print the distance matrix
    print_distance_matrix(num_coords, distance_matrix);

    // Free allocated memory
    for (int i = 0; i < num_coords; i++) {
        free(distance_matrix[i]);
    }
    free(distance_matrix);

    return 0;
}