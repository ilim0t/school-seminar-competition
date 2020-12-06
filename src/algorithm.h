#pragma once
#include "utils.h"

#include <stdbool.h>

void greedy_algorithm(const Param* const param,
                      const TSPdata* const tspdata,
                      Vdata* const vdata);

void greedy_depth_first_search(
    bool local_is_visited[],
    int minimum_cost[],
    int* local_tour,
    int local_depth,  // ここで探索するノードが何番目か、tour[local_depth]を順ぐりにまわす
    const Param* const param,
    const TSPdata* const tspdata,
    Vdata* const vdata);

void nearest_neighbor_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata);

void insertion_algorithm(const Param* const param,
                         const TSPdata* const tspdata,
                         Vdata* const vdata);

void two_approximation_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata);

int compute_min_upper_spanning_tree_kruskals(
    const TSPdata* const tspdata,
    bool adjacency_mat[tspdata->n][tspdata->n],
    int edge_dists[tspdata->n * (tspdata->n - 1) / 2][3]);

void prune_upper_spanning_tree(int n_prune_nodes,
                               const TSPdata* const tspdata,
                               bool adjacency_mat[tspdata->n][tspdata->n]);

void compute_edge_dists(int n,
                        int edge_dists[n * (n - 1) / 2][3],  // [dist, x, y]
                        const TSPdata* const tspdata);

void shaker_sort(int n, int array[n][3]);

void get_leaves(int root_node,
                int n,
                int* const node_in_tree_idx,  // node_in_tree[node_in_tree_idx]
                                              // = root_node と設定する
                int node_in_tree[n],
                bool const adjacency_mat[n][n]);

void shortcut_upper_spanning_tree(const TSPdata* const tspdata,
                                  bool adjacency_mat[tspdata->n][tspdata->n],
                                  int local_tour[tspdata->n]);

void stroke_depth_first_search(
    int root_node,
    int local_tour[],
    int* local_tour_idx,  // local_tour[*local_tour_idx] = root_node と設定する
    const TSPdata* const tspdata,
    bool adjacency_mat[tspdata->n][tspdata->n]);

void print_tour(const int* const local_tour, const TSPdata* const tspdata);
