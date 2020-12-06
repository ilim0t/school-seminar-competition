#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void two_approximation_algorithm(const Param* const param,
                                 const TSPdata* const tspdata,
                                 Vdata* const vdata) {
  bool adjacency_mat[tspdata->n][tspdata->n];
  for (int i = 0; i < tspdata->n; i++) {
    for (int j = 0; j < tspdata->n; j++) {
      adjacency_mat[i][j] = false;
    }
  }

  int edge_dists[tspdata->n * (tspdata->n - 1) / 2][3];
  compute_edge_dists(tspdata->n, edge_dists, tspdata);

  int n_node_spanning_tree = compute_min_upper_spanning_tree_kruskals(
      tspdata, adjacency_mat, edge_dists);
  prune_upper_spanning_tree(n_node_spanning_tree - tspdata->min_node_num,
                            tspdata, adjacency_mat);
#if DEBUG
  printf("[PRUNNED ADJACENCY MATRIX]\n");
  for (int i = 0; i < tspdata->n; i++) {
    for (int j = 0; j < tspdata->n; j++) {
      printf("%d ", adjacency_mat[i][j]);
    }
    printf("\n");
  }
  printf("\n");
#endif

  int min_cost = INT_MAX;

  while (cpu_time() - vdata->starttime < param->timelim) {
    int local_tour[tspdata->n];
    for (int i = 0; i < tspdata->n; i++) {
      local_tour[i] = -1;
    }
    shortcut_upper_spanning_tree(tspdata, adjacency_mat, local_tour);

    const int cost = compute_cost(tspdata, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
      for (int tour_idx = 0; tour_idx < tspdata->n; tour_idx++) {
        vdata->bestsol[tour_idx] = local_tour[tour_idx];
      }
    }
  }
}

void compute_edge_dists(int n,
                        int edge_dists[n * (n - 1) / 2][3],  // [dist, x, y]
                        const TSPdata* const tspdata) {
  int edge_dists_idx = 0;

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      edge_dists[edge_dists_idx][0] = dist(i, j);
      edge_dists[edge_dists_idx][1] = i;
      edge_dists[edge_dists_idx][2] = j;
      edge_dists_idx++;
    }
  }
}

int compute_min_upper_spanning_tree_kruskals(
    const TSPdata* const tspdata,
    bool adjacency_mat[tspdata->n][tspdata->n],
    int edge_dists[tspdata->n * (tspdata->n - 1) / 2][3]) {
  bool out_of_use_mat[tspdata->n][tspdata->n];
  for (int i = 0; i < tspdata->n; i++) {
    for (int j = 0; j < tspdata->n; j++) {
      out_of_use_mat[i][j] = false;
    }
  }

  shaker_sort(tspdata->n * (tspdata->n - 1) / 2, edge_dists);

  int n_node_spanning_tree = 1;
  for (int iter = 0; iter < tspdata->n * (tspdata->n - 1) / 2; iter++) {
    int node_1 = edge_dists[iter][1];
    int node_2 = edge_dists[iter][2];
    if (out_of_use_mat[node_1][node_2]) {
      continue;
    }

    int tree_1_idx = 0;
    int node_in_tree_1[tspdata->n];
    int tree_2_idx = 0;
    int node_in_tree_2[tspdata->n];

    get_leaves(node_1, tspdata->n, &tree_1_idx, node_in_tree_1, adjacency_mat);
    get_leaves(node_2, tspdata->n, &tree_2_idx, node_in_tree_2, adjacency_mat);

    for (int i = 0; i < tree_1_idx; i++) {
      for (int j = 0; j < tree_2_idx; j++) {
        out_of_use_mat[node_in_tree_1[i]][node_in_tree_2[j]] = true;
        out_of_use_mat[node_in_tree_2[j]][node_in_tree_1[i]] = true;
      }
    }

    adjacency_mat[node_1][node_2] = true;
    adjacency_mat[node_2][node_1] = true;

    n_node_spanning_tree++;

#if DEBUG
    printf("iter=%d (n_node_spanning_tree=%d)\n", iter, n_node_spanning_tree);
    for (int i = 0; i < tspdata->n; i++) {
      for (int j = 0; j < tspdata->n; j++) {
        printf("%d ", adjacency_mat[i][j]);
      }
      printf("\n");
    }
    printf("\n");
    for (int i = 0; i < tspdata->n; i++) {
      for (int j = 0; j < tspdata->n; j++) {
        printf("%d ", out_of_use_mat[i][j]);
      }
      printf("\n");
    }
    printf("\n\n");
#endif
    if (n_node_spanning_tree >= tspdata->min_node_num &&
        tree_1_idx + tree_2_idx == n_node_spanning_tree) {
      return n_node_spanning_tree;
    }
  }
}

void prune_upper_spanning_tree(int n_prune_nodes,
                               const TSPdata* const tspdata,
                               bool adjacency_mat[tspdata->n][tspdata->n]) {
  for (int prune_idx = 0; prune_idx < n_prune_nodes; prune_idx++) {
    int max_reduced_dist = 0;
    int prune_node;
    int prune_connected_node;

    for (int temp_prune_node = 0; temp_prune_node < tspdata->n;
         temp_prune_node++) {
      int n_branch = 0;
      int connected_node;

      for (int temp_connected_node = 0; temp_connected_node < tspdata->n;
           temp_connected_node++) {
        if (!adjacency_mat[temp_prune_node][temp_connected_node]) {
          continue;
        } else if (n_branch > 1) {
          break;
        }
        n_branch++;
        connected_node = temp_connected_node;
      }

      if (n_branch != 1) {
        continue;
      }

      const int reduced_dist = dist(temp_prune_node, connected_node);
      if (reduced_dist > max_reduced_dist) {
        max_reduced_dist = reduced_dist;
        prune_node = temp_prune_node;
        prune_connected_node = connected_node;
      }
    }

    adjacency_mat[prune_node][prune_connected_node] = false;
    adjacency_mat[prune_connected_node][prune_node] = false;
  }
}

void get_leaves(int root_node,
                int n,
                int* const node_in_tree_idx,  // node_in_tree[node_in_tree_idx]
                                              // = root_node と設定する
                int node_in_tree[n],
                bool const adjacency_mat[n][n]) {
  node_in_tree[*node_in_tree_idx] = root_node;
  (*node_in_tree_idx)++;

  for (int temp_next_node = 0; temp_next_node < n; temp_next_node++) {
    if (!adjacency_mat[root_node][temp_next_node]) {
      continue;
    }

    bool unreached = true;
    for (int already_counted_idx = 0; already_counted_idx < *node_in_tree_idx;
         already_counted_idx++) {
      if (temp_next_node == node_in_tree[already_counted_idx]) {
        unreached = false;
        break;
      }
    }
    if (unreached) {
      get_leaves(temp_next_node, n, node_in_tree_idx, node_in_tree,
                 adjacency_mat);
    }
  }
}

void shortcut_upper_spanning_tree(const TSPdata* const tspdata,
                                  bool adjacency_mat[tspdata->n][tspdata->n],
                                  int local_tour[tspdata->n]) {
  int node_in_tree_idx = 0;
  int node_in_tree[tspdata->min_node_num];

  for (int i = 0; i < tspdata->n; i++) {
    for (int j = 0; j < tspdata->n; j++) {
      if (adjacency_mat[i][j]) {
        node_in_tree[node_in_tree_idx] = i;
        node_in_tree_idx++;
        break;
      }
    }
  }

  const int first_node = node_in_tree[rand() % tspdata->min_node_num];
  int local_tour_idx = 0;
  stroke_depth_first_search(first_node, local_tour, &local_tour_idx, tspdata,
                            adjacency_mat);
}

void stroke_depth_first_search(
    int root_node,
    int local_tour[],
    int* local_tour_idx,  // local_tour[*local_tour_idx] = root_node と設定する
    const TSPdata* const tspdata,
    bool adjacency_mat[tspdata->n][tspdata->n]) {
  for (int tour_idx = 0; tour_idx < *local_tour_idx; tour_idx++) {
    if (root_node == local_tour[tour_idx]) {
      return;
    }
  }

  local_tour[*local_tour_idx] = root_node;
  (*local_tour_idx)++;
#ifdef DEBUG
  printf("[UPDATE!]\n");
  print_tour(local_tour, tspdata);
#endif

  for (int temp_next_node = 0; temp_next_node < tspdata->n; temp_next_node++) {
    if (adjacency_mat[root_node][temp_next_node]) {
      stroke_depth_first_search(temp_next_node, local_tour, local_tour_idx,
                                tspdata, adjacency_mat);
    }
  }
}
