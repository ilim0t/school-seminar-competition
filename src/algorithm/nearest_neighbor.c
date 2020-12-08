#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void nearest_neighbor_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata) {
  for (int i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
  }

  nearest_neighbor(tspdata->n, tspdata->min_node_num, param->timelim,
                   tspdata->x, tspdata->y, vdata->bestsol);
}

void nearest_neighbor(int n_nodes,
                      int n_min_nodes,
                      double timelim,
                      double x_coords[n_nodes],
                      double y_coords[n_nodes],
                      int best_tour[n_nodes]) {
  double starttime = cpu_time();
  int min_cost = INT_MAX;
  if (my_is_feasible(n_nodes, n_min_nodes, best_tour)) {
    min_cost = my_compute_tour_cost(n_nodes, x_coords, y_coords, best_tour);
  }

  while (cpu_time() - starttime < timelim) {
    bool is_visiteds[n_nodes];
    int local_tour[n_nodes];
    for (int i = 0; i < n_nodes; i++) {
      is_visiteds[i] = false;
      local_tour[i] = -1;
    }

    const int first_node = rand() % n_nodes;
    local_tour[0] = first_node;
    is_visiteds[first_node] = true;

    for (int depth = 1; depth < n_min_nodes; depth++) {
      int nearest_cost = INT_MAX;

      for (int temp_next_node = 0; temp_next_node < n_nodes; temp_next_node++) {
        if (is_visiteds[temp_next_node]) {
          continue;
        }

        const int cost =
            my_dist(x_coords, y_coords, local_tour[depth - 1], temp_next_node);
        if (cost < nearest_cost) {
          nearest_cost = cost;
          local_tour[depth] = temp_next_node;
        }
      }
      is_visiteds[local_tour[depth]] = true;
    }

    const int cost =
        my_compute_tour_cost(n_nodes, x_coords, y_coords, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
#ifdef DEBUG
      printf("\n[UPDATE] nearest_neighbor\n");
      print_tour(n_nodes, n_min_nodes, x_coords, y_coords, local_tour);
#endif

      for (int tour_idx = 0; tour_idx < n_nodes; tour_idx++) {
        best_tour[tour_idx] = local_tour[tour_idx];
      }
    }
  }
  // insertion(tspdata->n, tspdata->min_node_num, param->timelim * 0.1, tspdata->x,
  //           tspdata->y, vdata->bestsol);
}

void restricted_nearest_neighbor(int n_nodes,
                                 int n_min_nodes,
                                 double timelim,
                                 double x_coords[n_nodes],
                                 double y_coords[n_nodes],
                                 int best_tour[n_nodes],
                                 bool can_visit[n_nodes]) {
  double starttime = cpu_time();
  int min_cost = INT_MAX;
  if (my_is_feasible(n_nodes, n_min_nodes, best_tour)) {
    min_cost = my_compute_tour_cost(n_nodes, x_coords, y_coords, best_tour);
  }

  while (cpu_time() - starttime < timelim) {
    bool is_visiteds[n_nodes];
    int local_tour[n_nodes];
    for (int i = 0; i < n_nodes; i++) {
      is_visiteds[i] = !can_visit[i];
      local_tour[i] = -1;
    }

    int first_node;
    while (true) {
      first_node = rand() % n_nodes;
      if (!is_visiteds[first_node]) {
        break;
      }
    }
    local_tour[0] = first_node;
    is_visiteds[first_node] = true;

    for (int depth = 1; depth < n_min_nodes; depth++) {
      int nearest_cost = INT_MAX;

      for (int temp_next_node = 0; temp_next_node < n_nodes; temp_next_node++) {
        if (is_visiteds[temp_next_node]) {
          continue;
        }

        const int cost =
            my_dist(x_coords, y_coords, local_tour[depth - 1], temp_next_node);
        if (cost < nearest_cost) {
          nearest_cost = cost;
          local_tour[depth] = temp_next_node;
        }
      }
      is_visiteds[local_tour[depth]] = true;
    }

    const int cost =
        my_compute_tour_cost(n_nodes, x_coords, y_coords, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
#ifdef DEBUG
      printf("\n[UPDATE] nearest_neighbor\n");
      print_tour(n_nodes, n_min_nodes, x_coords, y_coords, local_tour);
#endif

      for (int tour_idx = 0; tour_idx < n_nodes; tour_idx++) {
        best_tour[tour_idx] = local_tour[tour_idx];
      }
    }
  }
}