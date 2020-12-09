#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void two_opt_replace(const int n_nodes,
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

  int count = 0;
  while (cpu_time() - starttime < timelim) {
    count++;
    // two-opt
    {
      const int a_tour_idx = rand() % n_min_nodes;
      const int b_tour_idx =
          (a_tour_idx + 2 + rand() % (n_min_nodes - 2)) % n_min_nodes;

      const int c_tour_idx = (a_tour_idx + 1) % n_min_nodes;
      const int d_tour_idx = (b_tour_idx + 1) % n_min_nodes;

      const int a_node = local_tour[a_tour_idx];
      const int b_node = local_tour[b_tour_idx];
      const int c_node = local_tour[c_tour_idx];
      const int d_node = local_tour[d_tour_idx];

      const int cost_diff = weighted_adjacency_mat[a_node][b_node] +
                            weighted_adjacency_mat[c_node][d_node] -
                            weighted_adjacency_mat[a_node][c_node] -
                            weighted_adjacency_mat[b_node][d_node];

      const double energy_diff = (double)cost_diff / min_cost * n_min_nodes;
      const double temp = pow(0.05, (cpu_time() - starttime) / timelim);
      const double acceptance_prob = exp(-energy_diff / temp);

      // if (cost_diff <= 0 || (double)rand() / RAND_MAX < acceptance_prob) {
      if (cost_diff <= 0) {
#if DEBUG > 1
        if (cost_diff > 0) {
          printf(
              "[Acceptance] two_opt, cost_diff: %d,\ttemp: %f,\tprob: "
              "%f,\tnow_cost: %d\n",
              cost_diff, temp, acceptance_prob,
              my_compute_tour_cost_mat(n_nodes, weighted_adjacency_mat,
                                       local_tour));
        }
#endif

        int new_tour[n_min_nodes];
        int new_tour_idx = 0;
#if DEBUG
        for (int i = 0; i < n_min_nodes; i++) {
          new_tour[i] = -1;
        }
#endif

        for (int tour_idx = d_tour_idx; tour_idx != a_tour_idx;
             tour_idx = (tour_idx + 1) % n_min_nodes) {
          new_tour[new_tour_idx] = local_tour[tour_idx];
          new_tour_idx++;
        }
        new_tour[new_tour_idx] = a_node;
        new_tour_idx++;

        for (int tour_idx = b_tour_idx; tour_idx != c_tour_idx;
             tour_idx = tour_idx == 0 ? n_min_nodes - 1 : (tour_idx - 1)) {
          new_tour[new_tour_idx] = local_tour[tour_idx];
          new_tour_idx++;
        }
        new_tour[new_tour_idx] = c_node;
        // new_tour_idx++;

        for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
          local_tour[tour_idx] = new_tour[tour_idx];
        }

        const int cost = my_compute_tour_cost_mat(
            n_nodes, weighted_adjacency_mat, local_tour);
        if (cost < min_cost) {
          min_cost = cost;
          for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
            best_tour[tour_idx] = local_tour[tour_idx];
          }

#if DEBUG
          printf("\n[UPDATE] two_opt(cost_diff=%d)\n", cost_diff);
          print_tour_mat(n_nodes, n_min_nodes, weighted_adjacency_mat,
                         local_tour);
          printf("\n");
#endif
        }
      }
    }

    // replace
    if ((double)rand() / RAND_MAX < 0.1) {
      int delete_node_idx_in_tour;
      int insert_idx;

      while (true) {
        delete_node_idx_in_tour = rand() % n_min_nodes;
        insert_idx = rand() % n_min_nodes;

        if ((delete_node_idx_in_tour + 1) % n_min_nodes != insert_idx) {
          break;
        }
      }

      int pre_insert;
      if (insert_idx == 0) {
        pre_insert = local_tour[n_min_nodes - 1];
      } else {
        pre_insert = local_tour[insert_idx - 1];
      }

      int following_insert;
      if (delete_node_idx_in_tour == insert_idx) {
        following_insert = local_tour[(insert_idx + 1) % n_min_nodes];
      } else {
        following_insert = local_tour[insert_idx];
      }

      int insert_node;
      int added_cost;
      if (delete_node_idx_in_tour == insert_idx) {
        while (true) {
          insert_node = rand() % n_nodes;
          if (!is_visiteds[insert_node]) {
            added_cost = weighted_adjacency_mat[pre_insert][insert_node] +
                         weighted_adjacency_mat[insert_node][following_insert];
            break;
          }
        }
      } else {
        int min_insert_cost = INT_MAX;
        for (int temp_insert_node = 0; temp_insert_node < n_nodes;
             temp_insert_node++) {
          if (is_visiteds[temp_insert_node]) {
            continue;
          }

          const int cost =
              weighted_adjacency_mat[pre_insert][temp_insert_node] +
              weighted_adjacency_mat[temp_insert_node][following_insert];
          if (cost < min_insert_cost) {
            min_insert_cost = cost;
            insert_node = temp_insert_node;
          }
        }
        added_cost = min_insert_cost -
                     weighted_adjacency_mat[pre_insert][following_insert];
      }

      const int pre_delete = local_tour[delete_node_idx_in_tour == 0
                                            ? n_min_nodes - 1
                                            : delete_node_idx_in_tour - 1];
      const int delete_node = local_tour[delete_node_idx_in_tour];
      const int following_delete =
          local_tour[(delete_node_idx_in_tour + 1) % n_min_nodes];

      int removed_cost;
      if (insert_idx == delete_node_idx_in_tour) {
        removed_cost = weighted_adjacency_mat[pre_delete][delete_node] +
                       weighted_adjacency_mat[delete_node][following_delete];
      } else {
        removed_cost = weighted_adjacency_mat[pre_delete][delete_node] +
                       weighted_adjacency_mat[delete_node][following_delete] -
                       weighted_adjacency_mat[pre_delete][following_delete];
      }
      const int cost_diff = added_cost - removed_cost;

      const double energy_diff = (double)cost_diff / min_cost * n_min_nodes;
      const double temp = pow(0.05, (cpu_time() - starttime) / timelim);
      const double acceptance_prob = exp(-energy_diff / temp);

      // if (cost_diff <= 0 || (double)rand() / RAND_MAX < acceptance_prob) {
      if (cost_diff <= 0) {
#if DEBUG > 1
        if (cost_diff > 0) {
          printf(
              "[Acceptance] replace,\tcost_diff: %d,\ttemp: %f,\tprob: "
              "%f,\tnow_cost: %d\n",
              cost_diff, temp, acceptance_prob,
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
            insert_idx = -1;
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

        const int cost = my_compute_tour_cost_mat(
            n_nodes, weighted_adjacency_mat, local_tour);
        if (cost < min_cost) {
          min_cost = cost;
          for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
            best_tour[tour_idx] = local_tour[tour_idx];
          }

#if DEBUG
          printf("\n[UPDATE] replace(cost_diff=%d)\n", cost_diff);
          print_tour_mat(n_nodes, n_min_nodes, weighted_adjacency_mat,
                         local_tour);
          printf("\n");
#endif
        }
      }
    }
  }
  printf("count: %d\n", count);
}
