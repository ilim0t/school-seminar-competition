#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>

void greedy_algorithm(const Param* const param,
                      const TSPdata* const tspdata,
                      Vdata* const vdata) {
  bool is_visiteds[tspdata->n];
  int local_tour[tspdata->n];

  for (int i = 0; i < tspdata->n; i++) {
    is_visiteds[i] = false;
  }

  int** weighted_adjacency_mat = int_d2array(tspdata->n, tspdata->n);
  compute_weighted_adjacency_mat(tspdata->n, tspdata->x, tspdata->y,
                                 weighted_adjacency_mat);

  int min_cost = INT_MAX;
  greedy_depth_first_search(
      tspdata->n, tspdata->min_node_num, 0, &min_cost, is_visiteds, local_tour,
      vdata->bestsol, vdata->starttime, param->timelim, weighted_adjacency_mat);
}

void greedy_depth_first_search(
    const int n_nodes,
    const int n_min_nodes,
    const int
        depth,  // ここで探索するノードが何番目か、tour[local_depth]を順ぐりにまわす
    int* min_cost,
    bool is_visiteds[],
    int local_tour[],
    int best_tour[],
    const int starttime,
    const int timelim,
    int** weighted_adjacency_mat) {
  if (cpu_time() - starttime > timelim) {
    return;
  }

  for (int temp_next_node = 0; temp_next_node < n_nodes; temp_next_node++) {
    if (is_visiteds[temp_next_node]) {
      continue;
    }

    local_tour[depth] = temp_next_node;

    if (depth + 1 == n_min_nodes) {
      if (n_min_nodes < n_nodes) {
        local_tour[n_min_nodes] = -1;
      }
      const int cost =
          my_compute_tour_cost_mat(n_nodes, weighted_adjacency_mat, local_tour);

      if (cost < *min_cost) {
#ifdef DEBUG
        printf("[UPDATE!]\n");
        print_tour_mat(n_nodes, n_min_nodes, weighted_adjacency_mat,
                       local_tour);
#endif
        *min_cost = cost;
        for (int tour_idx = 0; tour_idx <= n_min_nodes; tour_idx++) {
          best_tour[tour_idx] = local_tour[tour_idx];
        }
      }
      continue;
    }

    is_visiteds[temp_next_node] = true;
    greedy_depth_first_search(n_nodes, n_min_nodes, depth + 1, min_cost,
                              is_visiteds, local_tour, best_tour, starttime,
                              timelim, weighted_adjacency_mat);
    is_visiteds[temp_next_node] = false;
  }
}
