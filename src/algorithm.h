#pragma once
#include "utils.h"

#include <stdbool.h>

void my_algorithm(Param* param, TSPdata* tspdata, Vdata* vdata);

void DFS(bool* local_is_visited,
         int* minimum_cost,
         int* local_tour,
         int local_depth,  // ここで探索するノードが何番目か、
                           // tour[local_depth]を順ぐりにまわす
         Param* param,
         const TSPdata* const tspdata,
         Vdata* const vdata);