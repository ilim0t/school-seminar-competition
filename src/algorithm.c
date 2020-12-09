#include "algorithm.h"
#include "cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

#define PRINT_TOUR_LIM 20

void my_algorithm(const Param* const param,
                  const TSPdata* const tspdata,
                  Vdata* const vdata) {
  for (int i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
  }

  int** weighted_adjacency_mat = int_d2array(tspdata->n, tspdata->n);
  compute_weighted_adjacency_mat(tspdata->n, tspdata->x, tspdata->y,
                                 weighted_adjacency_mat);

  nearest_neighbor(tspdata->n, tspdata->min_node_num,
                   (param->timelim - cpu_time() + vdata->starttime) * 0.1,
                   weighted_adjacency_mat, vdata->bestsol);

  two_opt_replace(tspdata->n, tspdata->min_node_num,
                  param->timelim - cpu_time() + vdata->starttime,
                  weighted_adjacency_mat, vdata->bestsol);

  // const double iter_tim_lim = (param->timelim - cpu_time() +
  // vdata->starttime) /
  //                             ((double)param->timelim / 4);

  // while (cpu_time() - vdata->starttime < param->timelim) {
  //   if (tspdata->n > tspdata->min_node_num) {
  //     replace(tspdata->n, tspdata->min_node_num,
  //             fmin(iter_tim_lim * 0.5,
  //                  param->timelim - cpu_time() + vdata->starttime),
  //             weighted_adjacency_mat, vdata->bestsol);
  //   }
  //   two_opt(tspdata->n, tspdata->min_node_num,
  //           fmin(iter_tim_lim * 0.5,
  //                param->timelim - cpu_time() + vdata->starttime),
  //           weighted_adjacency_mat, vdata->bestsol);
  // }
};

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

int my_compute_tour_cost_mat(int n,
                             int** weighted_adjacency_mat,
                             const int* const tour) {
  int tour_idx, cost = 0;

  for (tour_idx = 0; tour_idx < n - 1; tour_idx++) {
    if (tour[tour_idx + 1] < 0) {
      break;
    }
    cost += weighted_adjacency_mat[tour[tour_idx]][tour[tour_idx + 1]];
  }
  cost += weighted_adjacency_mat[tour[tour_idx]][tour[0]];
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
  printf("time: %f\n", cpu_time());
  printf("is_feasible: %d\n", my_is_feasible(n, min_node_num, tour));
}

void print_tour_mat(int n,
                    int min_node_num,
                    int** weighted_adjacency_mat,
                    const int* const tour) {
  printf("cost: %d\n",
         my_compute_tour_cost_mat(n, weighted_adjacency_mat, tour));

  printf("tour: ");
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
  printf("time: %f\n", cpu_time());
  printf("is_feasible: %d\n", my_is_feasible(n, min_node_num, tour));
}

int** int_d2array(const int row, const int column) {
  int** const array = malloc(sizeof(int*) * row);
  if (array == NULL) {
    printf("メモリが確保できません\n");
    exit(1);
  }

  for (int i = 0; i < row; i++) {
    array[i] = malloc(sizeof(int) * column);
    if (array[i] == NULL) {
      printf("d2array[0]メモリが確保できません\n");
      exit(1);
    }
  }
  return array;
}

void int_d2free(int** array, int row) {
  for (int i = 0; i < row; i++) {
    free(array[i]);
  }
  free(array);
}

bool** bool_d2array(int row, int column) {
  bool** array = malloc(sizeof(bool*) * row);
  if (array == NULL) {
    printf("メモリが確保できません\n");
    exit(1);
  }

  for (int i = 0; i < row; i++) {
    array[i] = malloc(sizeof(bool) * column);
    if (array[i] == NULL) {
      printf("d2array[0]メモリが確保できません\n");
      exit(1);
    }
  }
  return array;
}

void bool_d2free(bool** array, int row) {
  for (int i = 0; i < row; i++) {
    free(array[i]);
  }
  free(array);
}
