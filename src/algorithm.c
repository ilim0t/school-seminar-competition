#include "algorithm.h"
#include "cpu_time.h"

#include <limits.h>

void my_algorithm(Param* param, TSPdata* tspdata, Vdata* vdata) {
  bool local_is_visited[tspdata->n];
  int local_tour[tspdata->min_node_num + 1];
  int local_depth = 0;

  for (int i = 0; i < tspdata->n; i++) {
    local_is_visited[i] = false;
  }

  for (int i = 0; i < tspdata->min_node_num; i++) {
    vdata->bestsol[i] = i;
  }
  vdata->bestsol[tspdata->min_node_num] = -1;
  int minimum_cost = INT_MAX;

  DFS(local_is_visited, &minimum_cost, local_tour, 0, param, tspdata, vdata);

  // while (cpu_time() - vdata->starttime < param->timelim) {
  // for (i = 0; i < tspdata->n; i++) {
  //   tmp = vdata->bestsol[i];
  //   vdata->bestsol[i] = vdata->bestsol[(i + 1) % (tspdata->n)];
  //   vdata->bestsol[(i + 1) % (tspdata->n)] = tmp;
  // }
  // }

  // param->timelim;
  // tspdata->min_node_num;
  // int n = tspdata->n;
  // tspdata->x[n - 1];
  // tspdata->y[n - 1];
  // vdata->starttime;
  // vdata->bestsol;
  // vdata->bestsol[0] = 2;
  // vdata->bestsol[1] = 1;
  // vdata->bestsol[2] = 0;
  // vdata->bestsol[3] = -1;
}
void DFS(bool* local_is_visited,
         int* minimum_cost,
         int* local_tour,
         int local_depth,  // ここで探索するノードが何番目か、
                           // tour[local_depth]を順ぐりにまわす
         Param* param,
         const TSPdata* const tspdata,
         Vdata* const vdata) {
  if (cpu_time() - vdata->starttime > param->timelim) {
    return;
  }

  for (int next_idx = 0; next_idx < tspdata->n; next_idx++) {
    if (local_is_visited[next_idx]) {
      continue;
    }

    local_tour[local_depth] = next_idx;

    if (local_depth + 1 == tspdata->min_node_num) {
      local_tour[tspdata->min_node_num] = -1;
      const int cost = compute_cost(tspdata, local_tour);
      // print_tour(local_tour, tspdata, vdata);

      if (cost < *minimum_cost) {
        printf("[UPDATE!]\n");
        print_tour(local_tour, tspdata, vdata);
        *minimum_cost = cost;
        for (int i = 0; i <= tspdata->min_node_num; i++) {
          vdata->bestsol[i] = local_tour[i];
        }
      }
      continue;
    }

    local_is_visited[next_idx] = true;
    DFS(local_is_visited, minimum_cost, local_tour, local_depth + 1, param,
        tspdata, vdata);
    local_is_visited[next_idx] = false;
  }
}

void print_tour(int* local_tour, TSPdata* tspdata, Vdata* vdata) {
  printf("%d:\t", compute_cost(tspdata, local_tour));
  for (int i = 0; i < tspdata->n; i++) {
    if (local_tour[i] < 0) {
      break;
    }
    printf("%d, ", local_tour[i]);
  }

  printf("\n");
}