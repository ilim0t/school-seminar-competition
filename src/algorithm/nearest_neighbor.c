#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

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
