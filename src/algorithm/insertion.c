#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void insertion_algorithm(const Param* const param,
                         const TSPdata* const tspdata,
                         Vdata* const vdata) {
  for (int i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
  }

  int** weighted_adjacency_mat = int_d2array(tspdata->n, tspdata->n);
  compute_weighted_adjacency_mat(tspdata->n, tspdata->x, tspdata->y,
                                 weighted_adjacency_mat);

  insertion(tspdata->n, tspdata->min_node_num, param->timelim,
            weighted_adjacency_mat, vdata->bestsol);
}

int insertion(const int n_nodes,
              const int n_min_nodes,
              const double timelim,
              int** weighted_adjacency_mat,
              int best_tour[n_nodes]) {
  double starttime = cpu_time();
  int min_cost = INT_MAX;
  if (my_is_feasible(n_nodes, n_min_nodes, best_tour)) {
    min_cost =
        my_compute_tour_cost_mat(n_nodes, weighted_adjacency_mat, best_tour);
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

    // int insert_location;
    for (int depth = 1; depth < n_min_nodes;
         depth++) {  // bestsol[depth]を設定する
      int min_cost_diff = INT_MAX;
      int next_node;
      int insert_loc;

      for (int temp_next_node = 0; temp_next_node < n_nodes; temp_next_node++) {
        if (is_visiteds[temp_next_node]) {
          continue;
        }

        for (int temp_insert_loc = 1; temp_insert_loc <= depth;
             temp_insert_loc++) {  // insert_location: sliceで使うものと一緒
          const int previous_node = local_tour[temp_insert_loc - 1];
          const int following_node = local_tour[temp_insert_loc % depth];

          const int cost_diff =
              weighted_adjacency_mat[previous_node][temp_next_node] +
              weighted_adjacency_mat[temp_next_node][following_node] -
              weighted_adjacency_mat[previous_node][following_node];
          if (cost_diff < min_cost_diff) {
            min_cost_diff = cost_diff;

            next_node = temp_next_node;
            insert_loc = temp_insert_loc;
          }
        }
      }

      // insert node
      for (int tour_idx = depth; tour_idx > insert_loc; tour_idx--) {
        local_tour[tour_idx] = local_tour[tour_idx - 1];
      }
      local_tour[insert_loc] = next_node;
      is_visiteds[next_node] = true;
    }

    const int cost =
        my_compute_tour_cost_mat(n_nodes, weighted_adjacency_mat, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
#ifdef DEBUG
      printf("\n[UPDATE] insertion\n");
      print_tour_mat(n_nodes, n_min_nodes, weighted_adjacency_mat, local_tour);
#endif

      for (int tour_idx = 0; tour_idx < n_nodes; tour_idx++) {
        best_tour[tour_idx] = local_tour[tour_idx];
      }
    }
  }
}
