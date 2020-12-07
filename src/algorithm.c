#include "algorithm.h"
#include "cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void swap(int a[3], int b[3]) {
  for (int i = 0; i < 3; i++) {
    int temp = a[i];
    a[i] = b[i];
    b[i] = temp;
  }
}

void shaker_sort(int n, int array[n][3]) {
  for (int i = 0; i < n; i++) {
    // head to tail
    for (int j = i + 1; j < n; j++) {
      if (array[j][0] < array[j - 1][0]) {
        swap(array[j], array[j - 1]);
      }
    }
    n--;

    // tail to head
    for (int k = n - 1; k > i; k--) {
      if (array[k][0] < array[k - 1][0]) {
        swap(array[k], array[k - 1]);
      }
    }
  }
}

int my_dist(double x_coords[], double y_coords[], int i, int j) {
  return sqrt((x_coords[i] - x_coords[j]) * (x_coords[i] - x_coords[j]) +
              (y_coords[i] - y_coords[j]) * (y_coords[i] - y_coords[j])) +
         0.5;
}

int my_compute_tour_cost(int n,
                         double x_coords[],
                         double y_coords[],
                         const int* const tour) {
  int tour_idx, cost = 0;

  for (tour_idx = 0; tour_idx < n - 1; tour_idx++) {
    if (tour[tour_idx + 1] < 0) {
      break;
    }
    cost += my_dist(x_coords, y_coords, tour[tour_idx], tour[tour_idx + 1]);
  }
  cost += my_dist(x_coords, y_coords, tour[tour_idx], tour[0]);
  return cost;
}

bool my_is_feasible(int n, int min_node_num, const int tour[n]) {
  int k, *visited, flag = 1, num_visited = 0;
  visited = (int*)malloc_e(n * sizeof(int));

  for (k = 0; k < n; k++)
    visited[k] = 0;
  for (k = 0; k < n; k++) {
    if (tour[k] < 0) {
      break;
    }
    if (tour[k] >= n) {
      flag = 0;
      break;
    }
    if (visited[tour[k]]) {
      flag = 0;
      break;
    } else {
      visited[tour[k]] = 1;
      num_visited++;
    }
  }
  if (num_visited < min_node_num)
    flag = 0;

  free(visited);
  /* if tour is feasible (resp., not feasible), then flag=1 (resp., 0) */
  return flag;
}

void print_tour(int n,
                int min_node_num,
                double x_coords[n],
                double y_coords[n],
                const int* const tour) {
  printf("cost: %d\n", my_compute_tour_cost(n, x_coords, y_coords, tour));

  printf("tour: ");
#define PRINT_TOUR_LIM 20
  int i;
  for (i = 0; i < n; i++) {
    if (tour[i] < 0) {
      break;
    }
    if (i == PRINT_TOUR_LIM) {
      printf(", ...");
      continue;
    } else if (i > PRINT_TOUR_LIM) {
      continue;
    }
    if (i != 0) {
      printf(", ");
    }
    printf("%d", tour[i]);
  }
  if (i > PRINT_TOUR_LIM) {
    printf(", %d", tour[i] < 0 ? tour[i - 1] : tour[i]);
  }

  printf("\nlength: %d\n", i);
  printf("is_feasible: %d\n", my_is_feasible(n, min_node_num, tour));
}