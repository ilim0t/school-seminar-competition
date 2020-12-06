#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>

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
