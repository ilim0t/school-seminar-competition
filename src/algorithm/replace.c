#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void replace_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata) {
  int i;
  for (i = 0; i < tspdata->min_node_num; i++) {
    vdata->bestsol[i] = i;
  }
  for (; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
  }

  int** weighted_adjacency_mat = int_d2array(tspdata->n, tspdata->n);
  compute_weighted_adjacency_mat(tspdata->n, tspdata->x, tspdata->y,
                                 weighted_adjacency_mat);

  replace(tspdata->n, tspdata->min_node_num, param->timelim,
          weighted_adjacency_mat, vdata->bestsol);
}

void replace(const int n_nodes,
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

  int local_tour[n_nodes];
  for (int i = 0; i < n_nodes; i++) {
    local_tour[i] = best_tour[i];
  }

  bool is_visiteds[n_nodes];
  for (int i; i < n_min_nodes; i++) {
    is_visiteds[i] = false;
  }
  for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
    is_visiteds[local_tour[tour_idx]] = true;
  }

  while (cpu_time() - starttime < timelim) {
    int delete_node_idx_in_tour;
    int insert_node;
    int insert_idx;

    while (true) {
      delete_node_idx_in_tour = rand() % n_min_nodes;
      insert_idx = rand() % n_min_nodes;
      insert_node = rand() % n_nodes;

      if ((delete_node_idx_in_tour + 1) % n_min_nodes == insert_idx) {
        continue;
      }
      if (insert_idx == delete_node_idx_in_tour) {
        if (!is_visiteds[insert_node]) {
          break;
        }
      } else {
        if (!is_visiteds[insert_node] ||
            insert_node == local_tour[delete_node_idx_in_tour]) {
          break;
        }
      }
    }

    const int delete_node = local_tour[delete_node_idx_in_tour];

    int reduced_cost;
    const int pre_insert =
        local_tour[insert_idx == 0 ? n_min_nodes - 1 : insert_idx - 1];
    const int following_insert =
        local_tour[insert_idx == delete_node_idx_in_tour
                       ? (insert_idx + 1) % n_min_nodes
                       : insert_idx];

    const int pre_delete =
        local_tour[delete_node_idx_in_tour == 0 ? n_min_nodes - 1
                                                : delete_node_idx_in_tour - 1];
    const int following_delete =
        local_tour[(delete_node_idx_in_tour + 1) % n_min_nodes];

    if (insert_idx == delete_node_idx_in_tour) {
      const int added_dist =
          weighted_adjacency_mat[pre_insert][insert_node] +
          weighted_adjacency_mat[insert_node][following_insert];
      const int delete_dist =
          weighted_adjacency_mat[pre_delete][delete_node] +
          weighted_adjacency_mat[delete_node][following_delete];

      reduced_cost = delete_dist - added_dist;
    } else {
      const int added_dist =
          weighted_adjacency_mat[pre_insert][insert_node] +
          weighted_adjacency_mat[insert_node][following_insert] -
          weighted_adjacency_mat[pre_insert][following_insert];
      const int delete_dist =
          weighted_adjacency_mat[pre_delete][delete_node] +
          weighted_adjacency_mat[delete_node][following_delete] -
          weighted_adjacency_mat[pre_delete][following_delete];

      reduced_cost = delete_dist - added_dist;
    }

    const double energy_diff =
        (double)reduced_cost / min_cost * n_min_nodes / 2;
    const double temp = pow(0.05, (cpu_time() - starttime) / timelim);
    const double acceptance_prob = exp(energy_diff / temp);

    if (reduced_cost > 0) {
    } else if ((double)rand() / RAND_MAX < acceptance_prob) {
      if (n_min_nodes > 1500) {
        continue;
      }
    } else {
      continue;
    }
#if DEBUG > 1
    if (reduced_cost < 0) {
      printf(
          "[Acceptance] replace,\treduced_cost: %d,\ttemp: %f,\tprob: "
          "%f,\tnow_cost: %d\n",
          reduced_cost, temp, acceptance_prob,
          my_compute_tour_cost_mat(n_nodes, weighted_adjacency_mat,
                                   local_tour));
    }
#endif
#if DEBUG

    int new_tour[n_nodes];

    for (int tour_idx = 0; tour_idx < n_nodes; tour_idx++) {
      new_tour[tour_idx] = -1;
    }
#else
    int new_tour[n_min_nodes];
#endif
    int new_tour_idx = 0;
    int old_tour_idx = 0;

    while (new_tour_idx < n_min_nodes) {
      if (old_tour_idx == insert_idx) {
        new_tour[new_tour_idx] = insert_node;
        insert_idx = -insert_idx - 1;
        new_tour_idx++;
      } else if (old_tour_idx == delete_node_idx_in_tour) {
        old_tour_idx++;
      } else {
        new_tour[new_tour_idx] = local_tour[old_tour_idx];
        new_tour_idx++;
        old_tour_idx++;
      }
    }

    for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
      local_tour[tour_idx] = new_tour[tour_idx];
    }

    is_visiteds[delete_node] = false;
    is_visiteds[insert_node] = true;

    const int cost =
        my_compute_tour_cost_mat(n_nodes, weighted_adjacency_mat, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
      for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
        best_tour[tour_idx] = local_tour[tour_idx];
      }
#if DEBUG
      printf("\n[UPDATE] replace(reduced_cost=%d)\n", reduced_cost);
      print_tour_mat(n_nodes, n_min_nodes, weighted_adjacency_mat, best_tour);
      printf("\n");
#endif
    }
  }
}
