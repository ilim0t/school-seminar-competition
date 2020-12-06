#include "algorithm.h"
#include "cpu_time.h"

#include <limits.h>
#include <math.h>
#define DEBUG 1

void my_algorithm(const Param* const param,
                  const TSPdata* const tspdata,
                  Vdata* const vdata) {
  bool local_is_visiteds[tspdata->n];
  int local_tour[tspdata->min_node_num + 1];
  int local_depth = 0;

  for (int i = 0; i < tspdata->n; i++) {
    local_is_visiteds[i] = false;
  }

  for (int i = 0; i < tspdata->min_node_num; i++) {
    vdata->bestsol[i] = i;
  }
  vdata->bestsol[tspdata->min_node_num] = -1;
  int minimum_cost = INT_MAX;

  DFS(local_is_visiteds, &minimum_cost, local_tour, 0, param, tspdata, vdata);

  // while (cpu_time() - vdata->starttime < param->timelim) {
  // for (i = 0; i < tspdata->n; i++) {
  //   tmp = vdata->bestsol[i];
  //   vdata->bestsol[i] = vdata->bestsol[(i + 1) % (tspdata->n)];
  //   vdata->bestsol[(i + 1) % (tspdata->n)] = tmp;
  // }
  // }
}

void DFS(
    bool* local_is_visiteds,
    int* minimum_cost,
    int* local_tour,
    int local_depth,  // ここで探索するノードが何番目か、tour[local_depth]を順ぐりにまわす
    const Param* const param,
    const TSPdata* const tspdata,
    Vdata* const vdata) {
  if (cpu_time() - vdata->starttime > param->timelim) {
    return;
  }

  for (int next_idx = 0; next_idx < tspdata->n; next_idx++) {
    if (local_is_visiteds[next_idx]) {
      continue;
    }

    local_tour[local_depth] = next_idx;

    if (local_depth + 1 == tspdata->min_node_num) {
      local_tour[tspdata->min_node_num] = -1;
      const int cost = compute_cost(tspdata, local_tour);

      if (cost < *minimum_cost) {
#ifdef DEBUG
        printf("[UPDATE!]\n");
        print_tour(local_tour, tspdata);
#endif
        *minimum_cost = cost;
        for (int i = 0; i <= tspdata->min_node_num; i++) {
          vdata->bestsol[i] = local_tour[i];
        }
      }
      continue;
    }

    local_is_visiteds[next_idx] = true;
    DFS(local_is_visiteds, minimum_cost, local_tour, local_depth + 1, param,
        tspdata, vdata);
    local_is_visiteds[next_idx] = false;
  }
}

void nearest_neighbor(const Param* const param,
                      const TSPdata* const tspdata,
                      Vdata* const vdata) {
  bool is_visiteds[tspdata->n];
  is_visiteds[0] = true;
  for (int i = 1; i < tspdata->n; i++) {
    is_visiteds[i] = false;
  }

  vdata->bestsol[0] = 0;
  for (int i = 1; i <= tspdata->min_node_num; i++) {
    vdata->bestsol[i] = -1;
  }

  for (int depth = 1; depth < tspdata->min_node_num; depth++) {
    int nearest_dist = INT_MAX;

    for (int next_node = 0; next_node < tspdata->n; next_node++) {
      if (is_visiteds[next_node]) {
        continue;
      }

      const int dist = dist(vdata->bestsol[depth - 1], next_node);
      if (dist < nearest_dist) {
        nearest_dist = dist;
        vdata->bestsol[depth] = next_node;
#ifdef DEBUG
        printf("[UPDATE!]\n");
        print_tour(vdata->bestsol, tspdata);
#endif
      }
    }
    is_visiteds[vdata->bestsol[depth]] = true;
  }
  vdata->bestsol[tspdata->min_node_num] = -1;
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