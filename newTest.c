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
int order[N];

void cheapestInsertion(int start_city) {
    int visited[N] = {0};
    int current_city = start_city;
    visited[start_city] = 1;

    int count = 1;
    int distance = 0;
    int path[N];
    path[0] = start_city;

    while (count < N) {
        int minCost = INT_MAX;
        int next_city = -1;

        #pragma omp parallel for reduction(min:minCost)
        for (int i = 0; i < N; i++) {
            if (!visited[i]) {
                #pragma omp critical
                {
                    for (int j = 0; j < count; j++) {
                        int cost = graph[path[j]][i] + graph[i][current_city] - graph[path[j]][current_city];
                        if (cost < minCost) {
                            minCost = cost;
                            next_city = i;
                        }
                    }
                }
            }
        }

        distance += minCost;
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

int main() {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        cheapestInsertion(i);
    }

    printf("Minimum distance: %d\n", min_distance);
    printf("Order of cities visited: ");
    for (int i = 0; i < N; i++) {
        printf("%d ", order[i]);
    }
    printf("\n");

    return 0;
}
