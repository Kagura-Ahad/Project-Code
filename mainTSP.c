#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdbool.h>

#define NUM_CITIES 4

typedef struct {
    char name[500];
    double x;
    double y;
} City;

// Define an enum for the choice of algorithm
typedef enum {
    BRUTE_FORCE,
    GREEDY
} Algorithm;

double distance(City a, City b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

void swap(City *a, City *b) {
    City temp = *a;
    *a = *b;
    *b = temp;
}

void write_log(char *filename, char *algorithm, City *path, int path_length, double total_distance)
{
    // Open the file in "w" mode to clear its contents
    FILE *file = fopen(filename, "a");

    // Check if the file was opened correctly
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return;
    }

    // Write to file the algorithm used
    fprintf(file, "Algorithm: %s\n", algorithm);

    //Write to file the total distance and path
    fprintf(file, "Total distance: %f\n", total_distance);
    fprintf(file, "Path: ");

    // Write to file the path
    for (int i = 0; i < path_length; i++) {
        fprintf(file, "%s ", path[i].name);
        if (i < path_length - 1)
        {
            fprintf(file, "-> ");
        }
        
    }
    fprintf(file, "\n\n");
    fclose(file);
}

void bruteForce(City *cities, int l, int r, City *best_path, double *best_distance) {
    
    if (l == r) {
        double total_distance = 0;
        for (int i = 0; i < r; i++) {
            total_distance += distance(cities[i], cities[i + 1]);
        }
        total_distance += distance(cities[r], cities[0]); // return to start
        if (total_distance < *best_distance) {
            *best_distance = total_distance;
            memcpy(best_path, cities, (r + 1) * sizeof(City));
            write_log("log.txt", "Brute Force", best_path, r + 1, total_distance);
        }
    } else {
        for (int i = l; i <= r; i++) {
            swap(&cities[l], &cities[i]);
            bruteForce(cities, l + 1, r, best_path, best_distance);
            swap(&cities[l], &cities[i]); // backtrack
        }
    }
}

void greedy(City *cities, City *best_path, double *best_distance) {
    bool visited[NUM_CITIES];
    for (int i = 0; i < NUM_CITIES; i++) {
        visited[i] = false;
    }

    // Start at the first city
    best_path[0] = cities[0];
    visited[0] = true;

    for (int i = 1; i < NUM_CITIES; i++) {
        double min_distance = DBL_MAX;
        int min_index = -1;

        // Find the nearest unvisited city
        for (int j = 0; j < NUM_CITIES; j++) {
            if (!visited[j]) {
                double dist = distance(best_path[i - 1], cities[j]);
                if (dist < min_distance) {
                    min_distance = dist;
                    min_index = j;
                }
            }
        }

        // Visit the nearest city
        best_path[i] = cities[min_index];
        visited[min_index] = true;
        *best_distance += min_distance;

        // Write the current path and total distance to the log file
        write_log("log.txt", "Greedy", best_path, i + 1, *best_distance);
    }

    // Return to the start
    *best_distance += distance(best_path[NUM_CITIES - 1], best_path[0]);
    write_log("log.txt", "Greedy", best_path, NUM_CITIES, *best_distance);
}

void tsp(City *cities, int num_cities, City *best_path, double *best_distance, Algorithm algorithm) {
    switch (algorithm) {
        case BRUTE_FORCE:
            bruteForce(cities, 0, num_cities - 1, best_path, best_distance);
            break;
        case GREEDY:
            greedy(cities, best_path, best_distance);
            break;
    }
}

int main() {
    City cities[NUM_CITIES];
    City best_path[NUM_CITIES];
    double best_distance = DBL_MAX;

    FILE *names_file = fopen("names.txt", "r");
    FILE *coords_file = fopen("cyl_XY_coords.txt", "r");
    
    char state[500];
    if (names_file == NULL || coords_file == NULL) {
        printf("Error opening file\n");
        return 1;
    }

    for (int i = 0; i < NUM_CITIES; i++) {
        if (fscanf(names_file, "%[^,], %[^\n]\n", cities[i].name, state) != 2) {
            printf("Error reading names.txt\n");
            return 1;
        }
        if (fscanf(coords_file, "%lf %lf\n", &cities[i].x, &cities[i].y) != 2) {
            printf("Error reading cyl_XY_coords.txt\n");
            return 1;
        }
    }

    fclose(names_file);
    fclose(coords_file);

    char filename[500] = "log.txt";

    // Open the file in "w" mode to clear its contents
    FILE *file = fopen(filename, "w");

    // Check if the file was opened correctly
    if (file == NULL) {
        printf("Error opening file %s\n", filename);
        return 1;
    }

    // Write to file the algorithm used
    fprintf(file, "Dataset contains %d cities\n\n", NUM_CITIES);
    fclose(file);

    clock_t start, end; // Declare start and end variables
    start = clock(); // Start the timer

    tsp(cities, NUM_CITIES, best_path, &best_distance, BRUTE_FORCE);

    end = clock(); // End the timer

    double time_taken = ((double)end - start) / CLOCKS_PER_SEC; // Calculate the time taken

    file = fopen("log.txt", "a"); // Open the file in "a" mode to append to it
    if (file == NULL) {
        printf("Error opening file log.txt\n");
        return 1;
    }
    // Write the time taken to the file
    fprintf(file, "Time taken: %f\n\n", time_taken); 

    //Append the unit of time taken to the file
    fprintf(file, "Unit of time taken: seconds\n\n");
    fclose(file); // Close the file

    printf("Best path: ");
    for (int i = 0; i < NUM_CITIES; i++) {
        printf("%s ", best_path[i].name);   
    }
    printf("\nTotal distance: %f\n", best_distance);

    return 0;
}
