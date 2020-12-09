#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void spanning_subtree_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata) {
  int** weighted_adjacency_mat = int_d2array(tspdata->n, tspdata->n);
  compute_weighted_adjacency_mat(tspdata->n, tspdata->x, tspdata->y,
                                 weighted_adjacency_mat);

  int** sorted_edge_idxs = int_d2array(tspdata->n * (tspdata->n - 1) / 2, 2);
  make_edge_idxs(tspdata->n, sorted_edge_idxs);
  sort_edge_idxs(tspdata->n, sorted_edge_idxs, weighted_adjacency_mat);

  bool** spanning_tree_adjacency_mat = bool_d2array(tspdata->n, tspdata->n);
  for (int i = 0; i < tspdata->n; i++) {
    for (int j = 0; j < tspdata->n; j++) {
      spanning_tree_adjacency_mat[i][j] = false;
    }
  }
  compute_min_spanning_tree_kruskals(tspdata->n, weighted_adjacency_mat,
                                     spanning_tree_adjacency_mat,
                                     sorted_edge_idxs);

  // prune_spanning_tree(tspdata->n, spanning_tree_adjacency_mat,
  //                     weighted_adjacency_mat);
  heuristics_prune_min_spanning_tree(
      tspdata->n, tspdata->n - tspdata->min_node_num,
      spanning_tree_adjacency_mat, weighted_adjacency_mat);

  for (int tour_idx = 0; tour_idx < tspdata->n; tour_idx++) {
    vdata->bestsol[tour_idx] = -1;
  }

  bool can_visit[tspdata->n];
  for (int i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
    can_visit[i] = false;
    for (int j = 0; j < tspdata->n; j++) {
      if (spanning_tree_adjacency_mat[i][j]) {
        can_visit[i] = true;
        break;
      }
    }
  }

  double remain_time = (param->timelim - cpu_time() + vdata->starttime);

  restricted_nearest_neighbor(tspdata->n, tspdata->min_node_num,
                              remain_time * 0.1, weighted_adjacency_mat,
                              vdata->bestsol, can_visit);

  // subtree2tour(tspdata->n, spanning_tree_adjacency_mat, vdata->bestsol);

  two_opt(tspdata->n, tspdata->min_node_num, remain_time * 0.3,
          weighted_adjacency_mat, vdata->bestsol);
  three_opt(tspdata->n, tspdata->min_node_num,
            (param->timelim - cpu_time() + vdata->starttime),
            weighted_adjacency_mat, vdata->bestsol);

  int_d2free(weighted_adjacency_mat, tspdata->n);
  int_d2free(sorted_edge_idxs, tspdata->n * (tspdata->n - 1) / 2);
  bool_d2free(spanning_tree_adjacency_mat, tspdata->n);
}

void compute_weighted_adjacency_mat(const int n_nodes,
                                    double x_coords[n_nodes],
                                    double y_coords[n_nodes],
                                    int** weighted_adjacency_mat) {
  for (int i = 0; i < n_nodes; i++) {
    weighted_adjacency_mat[i][i] = 0;

    for (int j = i + 1; j < n_nodes; j++) {
      const int cost = my_dist(x_coords, y_coords, i, j);
      weighted_adjacency_mat[i][j] = cost;
      weighted_adjacency_mat[j][i] = cost;
    }
  }
}

void make_edge_idxs(const int n_nodes, int** edge_idxs) {
  int edges_vec_idx = 0;

  for (int i = 0; i < n_nodes; i++) {
    for (int j = i + 1; j < n_nodes; j++) {
      edge_idxs[edges_vec_idx][0] = i;
      edge_idxs[edges_vec_idx][1] = j;
      edges_vec_idx++;
    }
  }
}

void sort_edge_idxs(const int n_nodes,
                    int** edge_idxs,
                    int** weighted_adjacency_mat  // const
) {
  int compare_edge(const void* a, const void* b) {
    const int a_dist = weighted_adjacency_mat[(*(int**)a)[0]][(*(int**)a)[1]];
    const int b_dist = weighted_adjacency_mat[(*(int**)b)[0]][(*(int**)b)[1]];
    return a_dist - b_dist;
  }
  qsort(edge_idxs, n_nodes * (n_nodes - 1) / 2, sizeof(int*), compare_edge);
}

void compute_min_spanning_tree_kruskals(const int n_nodes,
                                        int** weighted_adjacency_mat,
                                        bool** spanning_tree_adjacency_mat,
                                        int** sorted_edge_dists  // const
) {
  bool** out_of_use_mat = bool_d2array(n_nodes, n_nodes);

  for (int i = 0; i < n_nodes; i++) {
    for (int j = 0; j < n_nodes; j++) {
      out_of_use_mat[i][j] = false;
    }
  }

  int n_nodes_spanning_tree = 1;
  for (int iter = 0; iter < n_nodes * (n_nodes - 1) / 2; iter++) {
    int left_node = sorted_edge_dists[iter][0];
    int right_node = sorted_edge_dists[iter][1];

    if (out_of_use_mat[left_node][right_node]) {
      continue;
    }

    int left_tree_idx = 0;
    int nodes_in_left_tree[n_nodes];

    int right_tree_idx = 0;
    int nodes_in_right_tree[n_nodes];

    get_descendants(left_node, n_nodes, &left_tree_idx, nodes_in_left_tree,
                    spanning_tree_adjacency_mat);
    get_descendants(right_node, n_nodes, &right_tree_idx, nodes_in_right_tree,
                    spanning_tree_adjacency_mat);

    for (int i = 0; i < left_tree_idx; i++) {
      const int left_descendant = nodes_in_left_tree[i];

      for (int j = 0; j < right_tree_idx; j++) {
        const int right_descendant = nodes_in_right_tree[j];

        out_of_use_mat[left_descendant][right_descendant] = true;
        out_of_use_mat[right_descendant][left_descendant] = true;
      }
    }

    spanning_tree_adjacency_mat[left_node][right_node] = true;
    spanning_tree_adjacency_mat[right_node][left_node] = true;

    n_nodes_spanning_tree++;

#if DEBUG >= 2
    printf("iter=%d (n_nodes_spanning_tree=%d)\n", iter, n_nodes_spanning_tree);
    printf("[ADJACENCY_MAT]\n");
    for (int i = 0; i < n_nodes; i++) {
      for (int j = 0; j < n_nodes; j++) {
        printf("%d ", spanning_tree_adjacency_mat[i][j]);
      }
      printf("\n");
    }
    printf("\n[OUT_OF_USE_MAT]\n");
    for (int i = 0; i < n_nodes; i++) {
      for (int j = 0; j < n_nodes; j++) {
        printf("%d ", out_of_use_mat[i][j]);
      }
      printf("\n");
    }
    printf("\n\n");
#endif
    if (n_nodes_spanning_tree >= n_nodes) {
      bool_d2free(out_of_use_mat, n_nodes);
      return;
    }
  }
}

// void prune_spanning_tree(const int n_nodes,
//                          const int n_min_nodes,
//                          bool** spanning_tree_adjacency_mat,
//                          int** weighted_adjacency_mat) {
//   int n_nodes_in_sub_tree = n_nodes;
//   while (true) {

//   }
// }

void heuristics_prune_min_spanning_tree(const int n_nodes,
                                        const int n_prune_nodes,
                                        bool** spanning_tree_adjacency_mat,
                                        int** weighted_adjacency_mat  // const
) {
  int num_branches[n_nodes];
  for (int i = 0; i < n_nodes; i++) {
    num_branches[i] = 0;
    for (int j = 0; j < n_nodes; j++) {
      if (spanning_tree_adjacency_mat[i][j]) {
        num_branches[i]++;
      }
    }
  }

  for (int prune_idx = 0; prune_idx < n_prune_nodes; prune_idx++) {
    int max_reduced_dist = 0;
    int prune_node;
    int prune_connected_node;

    for (int temp_prune_node = 0; temp_prune_node < n_nodes;
         temp_prune_node++) {
      if (num_branches[temp_prune_node] != 1) {
        continue;
      }

      int connected_node;
      for (int temp_connected_node = 0; temp_connected_node < n_nodes;
           temp_connected_node++) {
        if (spanning_tree_adjacency_mat[temp_prune_node][temp_connected_node]) {
          connected_node = temp_connected_node;
        }
      }

      const int reduced_dist =
          weighted_adjacency_mat[temp_prune_node][connected_node];
      if (reduced_dist > max_reduced_dist) {
        max_reduced_dist = reduced_dist;
        prune_node = temp_prune_node;
        prune_connected_node = connected_node;
      }
    }

    num_branches[prune_node]--;
    num_branches[prune_connected_node]--;
    spanning_tree_adjacency_mat[prune_node][prune_connected_node] = false;
    spanning_tree_adjacency_mat[prune_connected_node][prune_node] = false;
  }
}

void get_descendants(const int parent_node,
                     const int n_nodes,
                     // nodes_in_tree[node_in_tree_idx] = parent_node と設定する
                     int* const node_in_tree_idx,
                     int nodes_in_tree[n_nodes],
                     bool** spanning_tree_adjacency_mat) {
  nodes_in_tree[*node_in_tree_idx] = parent_node;
  (*node_in_tree_idx)++;

  for (int temp_child_node = 0; temp_child_node < n_nodes; temp_child_node++) {
    if (!spanning_tree_adjacency_mat[parent_node][temp_child_node]) {
      continue;
    }

    // TODO: メモリに乗せる
    bool unreached = true;
    for (int already_counted_idx = 0; already_counted_idx < *node_in_tree_idx;
         already_counted_idx++) {
      if (temp_child_node == nodes_in_tree[already_counted_idx]) {
        unreached = false;
        break;
      }
    }
    if (unreached) {
      get_descendants(temp_child_node, n_nodes, node_in_tree_idx, nodes_in_tree,
                      spanning_tree_adjacency_mat);
    }
  }
}

void subtree2tour(const int n_nodes,
                  bool** spanning_tree_adjacency_mat,
                  int local_tour[n_nodes]) {
  bool included_sub_tree = false;

  int first_node;
  while (true) {
    first_node = rand() % n_nodes;

    for (int i = 0; i < n_nodes; i++) {
      if (spanning_tree_adjacency_mat[first_node][i]) {
        included_sub_tree = true;
        break;
      }
    }
    if (included_sub_tree) {
      break;
    }
  }

  int local_tour_idx = 0;
  stroke_depth_first_search(n_nodes, first_node, local_tour, &local_tour_idx,
                            spanning_tree_adjacency_mat);
}

void stroke_depth_first_search(
    const int n_nodes,
    const int root_node,
    int local_tour[n_nodes],
    int* local_tour_idx,  // local_tour[*local_tour_idx]
                          // = root_node と設定する
    bool** spanning_tree_adjacency_mat) {
  for (int tour_idx = 0; tour_idx < *local_tour_idx; tour_idx++) {
    if (root_node == local_tour[tour_idx]) {
      return;
    }
  }

  local_tour[*local_tour_idx] = root_node;
  (*local_tour_idx)++;

  for (int temp_next_node = 0; temp_next_node < n_nodes; temp_next_node++) {
    if (spanning_tree_adjacency_mat[root_node][temp_next_node]) {
      stroke_depth_first_search(n_nodes, temp_next_node, local_tour,
                                local_tour_idx, spanning_tree_adjacency_mat);
    }
  }
}
