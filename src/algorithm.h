#pragma once
#include "utils.h"

#include <stdbool.h>

#define DEBUG 1

void greedy_algorithm(const Param* const param,
                      const TSPdata* const tspdata,
                      Vdata* const vdata);

void greedy_depth_first_search(bool local_is_visited[],
                               int minimum_cost[],
                               int* local_tour,
                               int local_depth,
                               const Param* const param,
                               const TSPdata* const tspdata,
                               Vdata* const vdata);

void nearest_neighbor_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata);

void nearest_neighbor(int n_nodes,
                      int n_min_nodes,
                      double timelim,
                      double x_coords[n_nodes],
                      double y_coords[n_nodes],
                      int best_tour[n_nodes]);

void insertion_algorithm(const Param* const param,
                         const TSPdata* const tspdata,
                         Vdata* const vdata);

int insertion(int n_nodes,
              int n_min_nodes,
              double timelim,
              double x_coords[n_nodes],
              double y_coords[n_nodes],
              int best_tour[n_nodes]);

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
                        int edge_dists[n * (n - 1) / 2][3],
                        const TSPdata* const tspdata);

void get_leaves(int root_node,
                int n,
                int* const node_in_tree_idx,
                int node_in_tree[n],
                bool const adjacency_mat[n][n]);

void shortcut_upper_spanning_tree(const TSPdata* const tspdata,
                                  bool adjacency_mat[tspdata->n][tspdata->n],
                                  int local_tour[tspdata->n]);

void stroke_depth_first_search(int root_node,
                               int local_tour[],
                               int* local_tour_idx,
                               const TSPdata* const tspdata,
                               bool adjacency_mat[tspdata->n][tspdata->n]);

void two_opt_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata);

void two_opt(int n_nodes,
             int n_min_nodes,
             double timelim,
             double x_coords[n_nodes],
             double y_coords[n_nodes],
             int best_tour[n_nodes]);

void three_opt_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata);

void three_opt(int n_nodes,
             int n_min_nodes,
             double timelim,
             double x_coords[n_nodes],
             double y_coords[n_nodes],
             int best_tour[n_nodes]);

void shaker_sort(int n, int array[n][3]);

int my_dist(double x_coords[], double y_coords[], int i, int j);

int my_compute_tour_cost(int n,
                         double x_coords[],
                         double y_coords[],
                         const int* const tour);

void print_tour(int n,
                double x_coords[n],
                double y_coords[n],
                const int* const tour);
