#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void insertion_algorithm(const Param* const param,
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

    // int insert_location;
    for (int depth = 1; depth < tspdata->min_node_num;
         depth++) {  // bestsol[depth]を設定する
      int min_cost_diff = INT_MAX;
      int next_node;
      int insert_loc;

      for (int temp_next_node = 0; temp_next_node < tspdata->n;
           temp_next_node++) {
        if (is_visiteds[temp_next_node]) {
          continue;
        }

        for (int temp_insert_loc = 1; temp_insert_loc <= depth;
             temp_insert_loc++) {  // insert_location: sliceで使うものと一緒
          const int previous_node = local_tour[temp_insert_loc - 1];
          const int following_node = local_tour[temp_insert_loc % depth];
          const int cost_diff = dist(previous_node, temp_next_node) +
                                dist(temp_next_node, following_node) -
                                dist(previous_node, following_node);
          if (cost_diff < min_cost_diff) {
            min_cost_diff = cost_diff;

            next_node = temp_next_node;
            insert_loc = temp_insert_loc;
          }
        }
      }

      // insert node
      for (int tour_idx = depth; tour_idx > insert_loc; tour_idx--) {
        local_tour[tour_idx] = local_tour[tour_idx - 1];
      }
      local_tour[insert_loc] = next_node;
      is_visiteds[next_node] = true;

      // #ifdef DEBUG
      //       printf("[UPDATE!]\n");
      //       print_tour(tspdata->n, tspdata->x, tspdata->y, vdata->bestsol);
      // #endif
    }

    const int cost = compute_cost(tspdata, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
#ifdef DEBUG
      printf("[UPDATE!]\n");
      print_tour(tspdata->n, tspdata->x, tspdata->y, local_tour);
#endif

      for (int tour_idx = 0; tour_idx < tspdata->n; tour_idx++) {
        vdata->bestsol[tour_idx] = local_tour[tour_idx];
      }
    }
  }
}
