#pragma once
#include "utils.h"

#include <stdbool.h>

void my_algorithm(const Param* const param,
                  const TSPdata* const tspdata,
                  Vdata* const vdata);

void DFS(
    bool* local_is_visited,
    int* minimum_cost,
    int* local_tour,
    int local_depth,  // ここで探索するノードが何番目か、tour[local_depth]を順ぐりにまわす
    const Param* const param,
    const TSPdata* const tspdata,
    Vdata* const vdata);

void nearest_neighbor(const Param* const param,
                      const TSPdata* const tspdata,
                      Vdata* const vdata);

void print_tour(const int* const local_tour, const TSPdata* const tspdata);
