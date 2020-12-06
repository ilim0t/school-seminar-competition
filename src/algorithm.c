#include "algorithm.h"
#include "cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

// #define DEBUG 1

void greedy_algorithm(const Param* const param,
                      const TSPdata* const tspdata,
                      Vdata* const vdata) {
  bool is_visiteds[tspdata->n];
  int local_tour[tspdata->n];

  for (int i = 0; i < tspdata->n; i++) {
    is_visiteds[i] = false;
  }

  int min_cost = INT_MAX;
  greedy_depth_first_search(is_visiteds, local_tour, &min_cost, 0, param,
                            tspdata, vdata);
}

void greedy_depth_first_search(
    bool is_visiteds[],
    int local_tour[],

    int* min_cost,
    int depth,  // ここで探索するノードが何番目か、tour[local_depth]を順ぐりにまわす
    const Param* const param,
    const TSPdata* const tspdata,
    Vdata* const vdata) {
  if (cpu_time() - vdata->starttime > param->timelim) {
    return;
  }

  for (int temp_next_node = 0; temp_next_node < tspdata->n; temp_next_node++) {
    if (is_visiteds[temp_next_node]) {
      continue;
    }

    local_tour[depth] = temp_next_node;

    if (depth + 1 == tspdata->min_node_num) {
      if (tspdata->min_node_num < tspdata->n) {
        local_tour[tspdata->min_node_num] = -1;
      }
      const int cost = compute_cost(tspdata, local_tour);

      if (cost < *min_cost) {
#ifdef DEBUG
        printf("[UPDATE!]\n");
        print_tour(local_tour, tspdata);
#endif
        *min_cost = cost;
        for (int tour_idx = 0; tour_idx <= tspdata->min_node_num; tour_idx++) {
          vdata->bestsol[tour_idx] = local_tour[tour_idx];
        }
      }
      continue;
    }

    is_visiteds[temp_next_node] = true;
    greedy_depth_first_search(is_visiteds, min_cost, local_tour, depth + 1,
                              param, tspdata, vdata);
    is_visiteds[temp_next_node] = false;
  }
}

void nearest_neighbor_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata) {
  int min_cost = INT_MAX;

  while (cpu_time() - vdata->starttime < param->timelim) {
    bool is_visiteds[tspdata->n];
    int local_tour[tspdata->n];
    for (int i = 0; i < tspdata->n; i++) {
      is_visiteds[i] = false;
      local_tour[i] = -1;
    }

    const int first_node = rand() % tspdata->n;
    local_tour[0] = first_node;
    is_visiteds[first_node] = true;

    for (int depth = 1; depth < tspdata->min_node_num; depth++) {
      int nearest_cost = INT_MAX;

      for (int temp_next_node = 0; temp_next_node < tspdata->n;
           temp_next_node++) {
        if (is_visiteds[temp_next_node]) {
          continue;
        }

        const int cost = dist(local_tour[depth - 1], temp_next_node);
        if (cost < nearest_cost) {
          nearest_cost = cost;
          local_tour[depth] = temp_next_node;
          // #ifdef DEBUG
          //           printf("[UPDATE!]\n");
          //           print_tour(local_tour, tspdata);
          // #endif
        }
      }
      is_visiteds[local_tour[depth]] = true;
    }

    const int cost = compute_cost(tspdata, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
#ifdef DEBUG
      printf("[UPDATE!]\n");
      print_tour(local_tour, tspdata);
#endif

      for (int tour_idx = 0; tour_idx < tspdata->n; tour_idx++) {
        vdata->bestsol[tour_idx] = local_tour[tour_idx];
      }
    }
  }
}

void insertion_algorithm(const Param* const param,
                         const TSPdata* const tspdata,
                         Vdata* const vdata) {
  int min_cost = INT_MAX;

  while (cpu_time() - vdata->starttime < param->timelim) {
    bool is_visiteds[tspdata->n];
    int local_tour[tspdata->n];
    for (int i = 0; i < tspdata->n; i++) {
      is_visiteds[i] = false;
      local_tour[i] = -1;
    }

    const int first_node = rand() % tspdata->n;
    local_tour[0] = first_node;
    is_visiteds[first_node] = true;

    // int insert_location;
    for (int depth = 1; depth < tspdata->min_node_num;
         depth++) {  // bestsol[depth]を設定する
      int min_cost_diff = INT_MAX;
      int next_node;
      int insert_loc;

      for (int temp_next_node = 0; temp_next_node < tspdata->n;
           temp_next_node++) {
        if (is_visiteds[temp_next_node]) {
          continue;
        }

        for (int temp_insert_loc = 1; temp_insert_loc <= depth;
             temp_insert_loc++) {  // insert_location: sliceで使うものと一緒
          const int previous_node = local_tour[temp_insert_loc - 1];
          const int following_node = local_tour[temp_insert_loc % depth];
          const int cost_diff = dist(previous_node, temp_next_node) +
                                dist(temp_next_node, following_node) -
                                dist(previous_node, following_node);
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

      // #ifdef DEBUG
      //       printf("[UPDATE!]\n");
      //       print_tour(vdata->bestsol, tspdata);
      // #endif
    }

    const int cost = compute_cost(tspdata, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
#ifdef DEBUG
      printf("[UPDATE!]\n");
      print_tour(local_tour, tspdata);
#endif

      for (int tour_idx = 0; tour_idx < tspdata->n; tour_idx++) {
        vdata->bestsol[tour_idx] = local_tour[tour_idx];
      }
    }
  }
}

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

  int local_tour[tspdata->n];
  for (int i = 0; i < tspdata->n; i++) {
    local_tour[i] = -1;
  }
  shortcut_upper_spanning_tree(tspdata, adjacency_mat, local_tour);

  for (int tour_idx = 0; tour_idx < tspdata->n; tour_idx++) {
    vdata->bestsol[tour_idx] = local_tour[tour_idx];
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
  printf("%d\n", adjacency_mat[36][47]);
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

void print_tour(const int* const local_tour, const TSPdata* const tspdata) {
  printf("%d:\t", compute_cost(tspdata, local_tour));
  for (int i = 0; i < tspdata->n; i++) {
    if (local_tour[i] < 0) {
      break;
    }
    printf("%d, ", local_tour[i]);
  }

  printf("\n");
}