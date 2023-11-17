#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

#define N 4 // Number of cities

int graph[N][N] = {
    {0, 10, 15, 20},
    {10, 0, 35, 25},
    {15, 35, 0, 30},
    {20, 25, 30, 0}
};

int min_distance = INT_MAX;

void TSP(int start_city) {
    int visited[N] = {0};
    int current_city = start_city;
    visited[start_city] = 1;

    int count = 1;
    int distance = 0;

    while (count < N) {
        int nearest = INT_MAX;
        #pragma omp parallel for reduction(min:nearest)
        for (int i = 0; i < N; i++) {
            if (graph[current_city][i] < nearest && !visited[i]) {
                nearest = graph[current_city][i];
                start_city = i;
            }
        }

        distance += nearest;
        visited[start_city] = 1;
        current_city = start_city;
        count++;
    }

    distance += graph[current_city][start_city]; // Complete the cycle

    #pragma omp critical
    {
        if (distance < min_distance) {
            min_distance = distance;
        }
    }
}

int main() {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        TSP(i);
    }

    printf("Minimum distance: %d\n", min_distance);

    return 0;
}