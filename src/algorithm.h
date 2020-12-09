#pragma once
#include "utils.h"

#include <stdbool.h>

// #define DEBUG 1
// #define DEBUG 2

void my_algorithm(const Param* const param,
                  const TSPdata* const tspdata,
                  Vdata* const vdata);

// greedy.c
void greedy_algorithm(const Param* const param,
                      const TSPdata* const tspdata,
                      Vdata* const vdata);
void greedy_depth_first_search(const int n_nodes,
                               const int n_min_nodes,
                               const int depth,
                               int* min_cost,
                               bool is_visiteds[],
                               int local_tour[],
                               int best_tour[],
                               const int starttime,
                               const int timelim,
                               int** weighted_adjacency_mat);

// nearest_neighbor.c
void nearest_neighbor_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata);
void nearest_neighbor(const int n_nodes,
                      const int n_min_nodes,
                      const double timelim,
                      int** weighted_adjacency_mat,
                      int best_tour[n_nodes]);
void restricted_nearest_neighbor(const int n_nodes,
                                 const int n_min_nodes,
                                 const double timelim,
                                 int** weighted_adjacency_mat,
                                 int best_tour[n_nodes],
                                 bool can_visit[n_nodes]);

// insertion.c
void insertion_algorithm(const Param* const param,
                         const TSPdata* const tspdata,
                         Vdata* const vdata);
int insertion(const int n_nodes,
              const int n_min_nodes,
              const double timelim,
              int** weighted_adjacency_mat,
              int best_tour[n_nodes]);

// two_opt.c
void two_opt_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata);
void two_opt(const int n_nodes,
             const int n_min_nodes,
             const double timelim,
             int** weighted_adjacency_mat,
             int best_tour[n_nodes]);

// three_opt.c
void three_opt_algorithm(const Param* const param,
                         const TSPdata* const tspdata,
                         Vdata* const vdata);
void three_opt(const int n_nodes,
               const int n_min_nodes,
               const double timelim,
               int** weighted_adjacency_mat,
               int best_tour[n_nodes]);

// spanning_tree.c
void spanning_subtree_algorithm(const Param* const param,
                                const TSPdata* const tspdata,
                                Vdata* const vdata);
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
                                        int** weighted_adjacency_mat,
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
void replace(const int n_nodes,
             const int n_min_nodes,
             const double timelim,
             int** weighted_adjacency_matx_coords,
             int best_tour[n_nodes]);

// utils
int** int_d2array(const int row, const int column);
void int_d2free(int** array, const int row);
bool** bool_d2array(const int row, const int column);
void bool_d2free(bool** array, const int row);

int my_dist(double x_coords[], double y_coords[], const int i, const int j);
int my_compute_tour_cost(const int n,
                         double x_coords[],
                         double y_coords[],
                         const int* const tour);
int my_compute_tour_cost_mat(const int n,
                             int** weighted_adjacency_mat,
                             const int* const tour);
bool my_is_feasible(const int n, int min_node_num, const int tour[n]);

void print_tour(const int n,
                const int min_node_num,
                double x_coords[n],
                double y_coords[n],
                const int* const tour);
void print_tour_mat(const int n,
                    const int min_node_num,
                    int** weighted_adjacency_mat,
                    const int* const tour);
