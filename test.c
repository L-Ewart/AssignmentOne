#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include coordReader.c
#ifndef coordReader_H_
#define coordReader_H_

int readNumOfCoords();
double **readCoords ();
void *writeTourToFile();
#endif

int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);

// Struct to represent a point
struct Point {
    double x;
    double y;
};

// Function to calculate the Euclidean distance between two points
double calculate_distance(struct Point p1, struct Point p2) {
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

// Function to calculate the distance matrix for a given set of coordinates
void calculate_distance_matrix(struct Point coords[], int num_coords, double **distance_matrix) {
    for (int i = 0; i < num_coords; i++) {
        for (int j = 0; j < num_coords; j++) {
            distance_matrix[i][j] = calculate_distance(coords[i], coords[j]);
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
    // Call the function from Program 1 to generate the number of coordinates
    char input_filename[] = "9_coords.coord";
    int num_coords = readNumOfCoords(input_filename);
    printf("Number of coordinates: %d\n", num_coords);

    // Call the function from Program 1 to generate coordinates
    struct Point **coords = readCoords(input_filename, num_coords);

    // Dynamically allocate memory for the distance matrix
    double **distance_matrix = malloc(num_coords * sizeof(double *));
    for (int i = 0; i < num_coords; i++) {
        distance_matrix[i] = malloc(num_coords * sizeof(double));
    }

    // Calculate the distance matrix
    calculate_distance_matrix(coords, num_coords, distance_matrix);

    // Print the distance matrix
    print_distance_matrix(num_coords, distance_matrix);

    // Free allocated memory
    free(coords);
    for (int i = 0; i < num_coords; i++) {
        free(distance_matrix[i]);
    }
    free(distance_matrix);

    return 0;
}