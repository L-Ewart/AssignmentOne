#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#ifndef coordReader_H_
#define coordReader_H_

int readNumOfCoords();
double **readCoords ();
void *writeTourToFile();
#endif

int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);

int main(int argc, char *argv[]) {

    if (argc != 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1; // Return an error code
    }
    char *filename = argv[1];
    int readNumOfCoords(char *filename);
    double **readCoords(char *filename, int numOfCoords);
    // Access the file name from the command-line arguments
    //printf("Number of coordinates: %d\n", readNumOfCoords);
    return 0;
}