#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

#define N 9 // Number of cities

double graph[N][N] = {
    {0.000000,    676.464449,    258.496689,    583.277862,    221.195874,    353.277079,    453.909046,    446.711754,    807.890295,
    676.464449,    0.000000,    878.497983,    1112.471620,    712.726622,    959.411738,    256.758928,    689.716521,    671.811979,
    258.496689,   878.497983,    0.000000,    333.414398,    442.951092,    94.780396,    629.528098,    399.576868,    813.631469,
    583.277862,    1112.471620,    333.414398,    0.000000,    776.334790,    246.324495,    855.890183,    462.842661,    827.776548,
    221.195874,    712.726622,    442.951092,    776.334790,    0.000000,    532.613417,    549.540492,    661.263017,    993.881943,
    353.277079,    959.411738,    94.780396,    246.324495,    532.613417,    0.000000,    706.814888,    422.555358,    836.038955,
    453.909046,    256.758928,    629.528098,    855.890183,    549.540492,    706.814888,   0.000000,    445.308720,    542.271870,
    446.711754,    689.716521,    399.576868,    462.842661,    661.263017,    422.555358,    445.308720,    0.000000,    415.084260,
    807.890295,    671.811979,    813.631469,    827.776548,    993.881943,    836.038955,    542.271870,    415.084260,    0.000000}
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
