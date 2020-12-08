#pragma once
#include "utils.h"

#include <stdbool.h>

#define DEBUG 1
// #define DEBUG 2

// greedy.c
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

// nearest_neighbor.c
void nearest_neighbor_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata);
void nearest_neighbor(int n_nodes,
                      int n_min_nodes,
                      double timelim,
                      double x_coords[n_nodes],
                      double y_coords[n_nodes],
                      int best_tour[n_nodes]);
void restricted_nearest_neighbor(int n_nodes,
                                 int n_min_nodes,
                                 double timelim,
                                 double x_coords[n_nodes],
                                 double y_coords[n_nodes],
                                 int best_tour[n_nodes],
                                 bool can_visit[n_nodes]);

// insertion.c
void insertion_algorithm(const Param* const param,
                         const TSPdata* const tspdata,
                         Vdata* const vdata);
int insertion(int n_nodes,
              int n_min_nodes,
              double timelim,
              double x_coords[n_nodes],
              double y_coords[n_nodes],
              int best_tour[n_nodes]);

// 仮
void two_approximation_algorithm(const Param* const param,
                                 const TSPdata* const tspdata,
                                 Vdata* const vdata);

// two_opt.c
void two_opt_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata);
void two_opt(int n_nodes,
             int n_min_nodes,
             double timelim,
             double x_coords[n_nodes],
             double y_coords[n_nodes],
             int best_tour[n_nodes]);

// three_opt.c
void three_opt_algorithm(const Param* const param,
                         const TSPdata* const tspdata,
                         Vdata* const vdata);
void three_opt(int n_nodes,
               int n_min_nodes,
               double timelim,
               double x_coords[n_nodes],
               double y_coords[n_nodes],
               int best_tour[n_nodes]);

// spanning_tree.c
void compute_weighted_adjacency_mat(const int n_nodes,
                                    double x_coords[n_nodes],
                                    double y_coords[n_nodes],
                                    int** weighted_adjacency_mat);

void make_edge_idxs(const int n_nodes, int** edge_idxs);
void sort_edge_idxs(const int n_nodes,
                    int** edge_idxs,
                    int** weighted_adjacency_mat  // const
);
void sort_edge_idxs(const int n_nodes,
                    int** edge_idxs,
                    int** weighted_adjacency_mat  // const
);
void compute_min_spanning_tree_kruskals(const int n_nodes,
                                        double x_coords[n_nodes],
                                        double y_coords[n_nodes],
                                        bool** adjacency_mat,
                                        int** sorted_edge_dists  // const
);
void get_descendants(const int parent_node,
                     const int n_nodes,
                     int* const node_in_tree_idx,
                     int nodes_in_tree[n_nodes],
                     bool** adjacency_mat);

void heuristics_prune_min_spanning_tree(const int n_nodes,
                                        const int n_prune_nodes,
                                        bool** spanning_tree_adjacency_mat,
                                        int** weighted_adjacency_mat);
void subtree2tour(const int n_nodes,
                  bool** spanning_tree_adjacency_mat,
                  int local_tour[n_nodes]);
void stroke_depth_first_search(const int n_nodes,
                               const int root_node,
                               int local_tour[n_nodes],
                               int* local_tour_idx,
                               bool** spanning_tree_adjacency_mat);

// replace
void replace_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata);
void replace(int n_nodes,
             int n_min_nodes,
             double timelim,
             double x_coords[n_nodes],
             double y_coords[n_nodes],
             int best_tour[n_nodes]);

// utils
int** int_d2array(int row, int column);
void int_d2free(int** array, int row);
bool** bool_d2array(int row, int column);
void bool_d2free(bool** array, int row);

int my_dist(double x_coords[], double y_coords[], int i, int j);
int my_compute_tour_cost(int n,
                         double x_coords[],
                         double y_coords[],
                         const int* const tour);
bool my_is_feasible(int n, int min_node_num, const int tour[n]);

void print_tour(int n,
                int min_node_num,
                double x_coords[n],
                double y_coords[n],
                const int* const tour);

void spanning_subtree_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata);
