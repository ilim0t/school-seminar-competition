#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void two_opt_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata) {
  for (int i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
  }
  nearest_neighbor(tspdata->n, tspdata->min_node_num, param->timelim * 0.1,
                   tspdata->x, tspdata->y, vdata->bestsol);
  insertion(tspdata->n, tspdata->min_node_num, param->timelim * 0.1, tspdata->x,
            tspdata->y, vdata->bestsol);

#if DEBUG
  printf("\nnearest_neighbor\n");
  print_tour(tspdata->n, tspdata->x, tspdata->y, vdata->bestsol);
#endif

  while (cpu_time() - vdata->starttime < param->timelim) {
    const int a_tour_idx = rand() % tspdata->min_node_num;
    const int c_tour_idx = (a_tour_idx + 1) % tspdata->min_node_num;

    const int b_tour_idx =
        (a_tour_idx + 1 + rand() % (tspdata->min_node_num - 1)) %
        tspdata->min_node_num;
    const int d_tour_idx = (b_tour_idx + 1) % tspdata->min_node_num;

    const int a_node = vdata->bestsol[a_tour_idx];
    const int b_node = vdata->bestsol[b_tour_idx];
    const int c_node = vdata->bestsol[c_tour_idx];
    const int d_node = vdata->bestsol[d_tour_idx];

    // if (a < d) {

    // }
    const int reduces_cost = my_dist(tspdata->x, tspdata->y, a_node, c_node) +
                             my_dist(tspdata->x, tspdata->y, b_node, d_node) -
                             my_dist(tspdata->x, tspdata->y, a_node, b_node) -
                             my_dist(tspdata->x, tspdata->y, c_node, d_node);
    if (reduces_cost <= 0) {
      continue;
    }

    int new_tour[tspdata->min_node_num];
    int new_tour_idx = 0;
#if DEBUG
    for (int i = 0; i < tspdata->min_node_num; i++) {
      new_tour[i] = -1;
    }
#endif

    for (int tour_idx = b_tour_idx; tour_idx != c_tour_idx;
         tour_idx = tour_idx == 0 ? tspdata->min_node_num - 1
                                  : (tour_idx - 1)) {
      new_tour[new_tour_idx] = vdata->bestsol[tour_idx];
      new_tour_idx++;
    }
    new_tour[new_tour_idx] = vdata->bestsol[c_tour_idx];
    new_tour_idx++;

    for (int tour_idx = d_tour_idx; tour_idx != a_tour_idx;
         tour_idx = (tour_idx + 1) % tspdata->min_node_num) {
      new_tour[new_tour_idx] = vdata->bestsol[tour_idx];
      new_tour_idx++;
    }
    new_tour[new_tour_idx] = vdata->bestsol[a_tour_idx];
    new_tour_idx++;

    for (int tour_idx = 0; tour_idx < tspdata->min_node_num; tour_idx++) {
      vdata->bestsol[tour_idx] = new_tour[tour_idx];
    }
#if DEBUG
    printf("\n[UPDATE] (reduces_cost=%d)\n", reduces_cost);
    print_tour(tspdata->n, tspdata->x, tspdata->y, vdata->bestsol);
    printf("\n");
    printf(" ");
#endif
  }
}