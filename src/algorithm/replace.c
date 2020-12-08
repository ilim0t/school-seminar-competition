#include "../algorithm.h"
#include "../cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void replace_algorithm(const Param* const param,
                       const TSPdata* const tspdata,
                       Vdata* const vdata) {
  int i;
  for (i = 0; i < tspdata->min_node_num; i++) {
    vdata->bestsol[i] = i;
  }
  for (; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
  }

  replace(tspdata->n, tspdata->min_node_num, param->timelim, tspdata->x,
          tspdata->y, vdata->bestsol);
}

void replace(const int n_nodes,
             const int n_min_nodes,
             const double timelim,
             double x_coords[n_nodes],
             double y_coords[n_nodes],
             int best_tour[n_nodes]) {
  double starttime = cpu_time();
  int min_cost = INT_MAX;
  if (my_is_feasible(n_nodes, n_min_nodes, best_tour)) {
    min_cost = my_compute_tour_cost(n_nodes, x_coords, y_coords, best_tour);
  }

  int local_tour[n_nodes];
  for (int i = 0; i < n_nodes; i++) {
    local_tour[i] = best_tour[i];
  }

  while (cpu_time() - starttime < timelim) {
    int delete_node_idx_in_tour;
    int insert_node;
    int insert_idx;

    while (true) {
      delete_node_idx_in_tour = rand() % n_min_nodes;
      insert_node = rand() % n_nodes;
      insert_idx = rand() % n_min_nodes;

      if ((delete_node_idx_in_tour + 1) % n_min_nodes == insert_idx) {
        continue;
      }
      bool is_visited = false;
      for (int i = 0; i < n_min_nodes; ++i) {
        if (local_tour[i] == insert_node) {
          is_visited = true;
          break;
        }
      }
      if (!is_visited) {
        break;
      }
    }

    const int delete_node = local_tour[delete_node_idx_in_tour];

    int reduced_dist;
    const int pre_insert =
        local_tour[insert_idx == 0 ? n_min_nodes - 1 : insert_idx - 1];
    const int following_insert =
        local_tour[insert_idx == delete_node_idx_in_tour
                       ? (insert_idx + 1) % n_min_nodes
                       : insert_idx];

    const int pre_delete =
        local_tour[delete_node_idx_in_tour == 0 ? n_min_nodes - 1
                                                : delete_node_idx_in_tour - 1];
    const int following_delete =
        local_tour[(delete_node_idx_in_tour + 1) % n_min_nodes];

    if (insert_idx == delete_node_idx_in_tour) {
      const int added_dist =
          my_dist(x_coords, y_coords, pre_insert, insert_node) +
          my_dist(x_coords, y_coords, insert_node, following_insert);
      const int delete_dist =
          my_dist(x_coords, y_coords, pre_delete, delete_node) +
          my_dist(x_coords, y_coords, delete_node, following_delete);

      reduced_dist = delete_dist - added_dist;
    } else {
      const int added_dist =
          my_dist(x_coords, y_coords, pre_insert, insert_node) +
          my_dist(x_coords, y_coords, insert_node, following_insert) -
          my_dist(x_coords, y_coords, pre_insert, following_insert);
      const int delete_dist =
          my_dist(x_coords, y_coords, pre_delete, delete_node) +
          my_dist(x_coords, y_coords, delete_node, following_delete) -
          my_dist(x_coords, y_coords, pre_delete, following_delete);

      reduced_dist = delete_dist - added_dist;
    }

    if (reduced_dist < 0) {
      continue;
    }
#ifdef DEBUG
    int new_tour[n_nodes];

    for (int tour_idx = 0; tour_idx < n_nodes; tour_idx++) {
      new_tour[tour_idx] = -1;
    }
#else
    int new_tour[n_min_nodes];
#endif
    int new_tour_idx = 0;
    int old_tour_idx = 0;

    while (new_tour_idx < n_min_nodes) {
      if (old_tour_idx == insert_idx) {
        new_tour[new_tour_idx] = insert_node;
        insert_idx = -1;
        new_tour_idx++;
      } else if (old_tour_idx == delete_node_idx_in_tour) {
        old_tour_idx++;
      } else {
        new_tour[new_tour_idx] = local_tour[old_tour_idx];
        new_tour_idx++;
        old_tour_idx++;
      }
    }

    for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
      local_tour[tour_idx] = new_tour[tour_idx];
    }

    for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
      best_tour[tour_idx] = local_tour[tour_idx];
    }
    // const int cost = my_compute_tour_cost(n_nodes, x_coords, y_coords, )
    // const int cost = my_compyute()
    // if (min_cost )

#ifdef DEBUG
    printf("\n[UPDATE] replace(reduced=%d)\n", reduced_dist);
    print_tour(n_nodes, n_min_nodes, x_coords, y_coords, best_tour);
#endif
  }
}