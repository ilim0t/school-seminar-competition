#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void two_opt_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata) {
  for (int i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
  }

  int** weighted_adjacency_mat = int_d2array(tspdata->n, tspdata->n);
  compute_weighted_adjacency_mat(tspdata->n, tspdata->x, tspdata->y,
                                 weighted_adjacency_mat);

  nearest_neighbor(tspdata->n, tspdata->min_node_num, param->timelim * 0.2,
                   weighted_adjacency_mat, vdata->bestsol);

  two_opt(tspdata->n, tspdata->min_node_num,
          param->timelim - cpu_time() - vdata->starttime,
          weighted_adjacency_mat, vdata->bestsol);
}

void two_opt(const int n_nodes,
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

  while (cpu_time() - starttime < timelim) {
    const int a_tour_idx = rand() % n_min_nodes;
    const int c_tour_idx = (a_tour_idx + 1) % n_min_nodes;

    const int b_tour_idx =
        (a_tour_idx + 1 + rand() % (n_min_nodes - 1)) % n_min_nodes;
    const int d_tour_idx = (b_tour_idx + 1) % n_min_nodes;

    const int a_node = best_tour[a_tour_idx];
    const int b_node = best_tour[b_tour_idx];
    const int c_node = best_tour[c_tour_idx];
    const int d_node = best_tour[d_tour_idx];

    const int reduced_cost = weighted_adjacency_mat[a_node][c_node] +
                             weighted_adjacency_mat[b_node][d_node] -
                             weighted_adjacency_mat[a_node][b_node] -
                             weighted_adjacency_mat[c_node][d_node];

    const double energy_diff =
        (double)reduced_cost / min_cost * n_min_nodes / 2;
    const double temp = pow(0.05, (cpu_time() - starttime) / timelim);
    const double acceptance_prob = exp(energy_diff / temp);

    if (reduced_cost > 0) {
    } else if ((double)rand() / RAND_MAX < acceptance_prob) {
      if (n_min_nodes > 1000) {
        continue;
      }
    } else {
      continue;
    }
#if DEBUG > 1
    if (reduced_cost < 0) {
      printf("[Acceptance] two_opt, reduced_cost: %d, temp: %f, prob: %f\n",
             reduced_cost, temp, acceptance_prob);
    }
#endif

    int new_tour[n_min_nodes];
    int new_tour_idx = 0;
#if DEBUG
    for (int i = 0; i < n_min_nodes; i++) {
      new_tour[i] = -1;
    }
#endif

    for (int tour_idx = b_tour_idx; tour_idx != c_tour_idx;
         tour_idx = tour_idx == 0 ? n_min_nodes - 1 : (tour_idx - 1)) {
      new_tour[new_tour_idx] = best_tour[tour_idx];
      new_tour_idx++;
    }
    new_tour[new_tour_idx] = c_node;
    new_tour_idx++;

    for (int tour_idx = d_tour_idx; tour_idx != a_tour_idx;
         tour_idx = (tour_idx + 1) % n_min_nodes) {
      new_tour[new_tour_idx] = best_tour[tour_idx];
      new_tour_idx++;
    }
    new_tour[new_tour_idx] = a_node;
    // new_tour_idx++;

    for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
      local_tour[tour_idx] = new_tour[tour_idx];
    }

    const int cost =
        my_compute_tour_cost_mat(n_nodes, weighted_adjacency_mat, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
      for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
        best_tour[tour_idx] = local_tour[tour_idx];
      }
#if DEBUG
      printf("\n[UPDATE] two_opt(reduced_cost=%d)\n", reduced_cost);
      print_tour_mat(n_nodes, n_min_nodes, weighted_adjacency_mat, best_tour);
      printf("\n");
#endif
    }
  }
}