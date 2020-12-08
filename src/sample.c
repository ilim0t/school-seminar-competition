/*****************************************************************************
  A template program for 2-dimensional euclidean symmetric TSP solver.
  Subroutines to read instance data and compute the objective value of a given
  tour (solution) are included. The one to output the computed tour in the
  TSPLIB format is also included.

  The URL of TSPLIB is:
       http://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/

  NOTE: Indices of nodes range from 0 to n-1 in this program,
        while it does from 1 to n in "README.eng" and the data files of
        instances and of tours.

  If you would like to use various parameters, it might be useful to modify
  the definition of struct "Param" and mimic the way the default value of
  "timelim" is given and how its value is input from the command line.
******************************************************************************/

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cpu_time.c"

/***** constants *************************************************************/
#define MAX_STR 1024

/***** macros ****************************************************************/
#define dist(k, l)                                   \
  ((int)(sqrt((tspdata->x[k] - tspdata->x[l]) *      \
                  (tspdata->x[k] - tspdata->x[l]) +  \
              (tspdata->y[k] - tspdata->y[l]) *      \
                  (tspdata->y[k] - tspdata->y[l])) + \
         0.5))

/***** default values of parameters ******************************************/
#define TIMELIM 300 /* the time limit for the algorithm in seconds */
#define GIVESOL 0   /* 1: input a solution; 0: do not give a solution */
#define OUTFORMAT                                    \
  2 /* 0: do not output the tour;                    \
       1: output the computed tour in TSPLIB format; \
       2: output the computed tour in TSP_VIEW format */
#define TOURFILE "result.tour"
/* the output file of computed tour */

typedef struct {
  int timelim;            /* the time limit for the algorithm in secs. */
  int givesol;            /* give a solution (1) or not (0) */
  int outformat;          /* 1: output the computed tour in TSPLIB format;
                             0: do not output it */
  char tourfile[MAX_STR]; /* the output file of computed tour */
  /* NEVER MODIFY THE ABOVE VARIABLES.  */
  /* You can add more components below. */

} Param; /* parameters */

typedef struct {
  char name[MAX_STR]; /* name of the instance */
  int n;              /* number of nodes */
  double* x;          /* x-coordinates of nodes */
  double* y;          /* y-coordinates of nodes */
  int min_node_num;   /* minimum number of nodes the solution contains */
} TSPdata;            /* data of TSP instance */

typedef struct {
  double timebrid;  /* the time before reading the instance data */
  double starttime; /* the time the search started */
  double endtime;   /* the time the search ended */
  int* bestsol;     /* the best solution found so far */
  /* NEVER MODIFY THE ABOVE FOUR VARIABLES. */
  /* You can add more components below. */

} Vdata; /* various data often necessary during the search */

/************************ declaration of functions ***************************/
FILE* open_file(char* fname, char* mode);
void* malloc_e(size_t size);

void copy_parameters(int argc, char* argv[], Param* param);
void prepare_memory(TSPdata* tspdata, Vdata* vdata);
void read_header(FILE* in, TSPdata* tspdata);
void read_tspfile(FILE* in, TSPdata* tspdata, Vdata* vdata);
void read_tourfile(FILE* in, TSPdata* tspdata, int* tour);
void output_tour(FILE* out, TSPdata* tspdata, int* tour);
void output_tour_for_tsp_view(FILE* out, TSPdata* tspdata, int* tour);
void recompute_obj(Param* param, TSPdata* tspdata, Vdata* vdata);
void my_algorithm(Param* param, TSPdata* tspdata, Vdata* vdata);
void nearest_neighbor(TSPdata* tspdata,
                      const double timelim,
                      int best_tour[tspdata->n]);
void replace(TSPdata* tspdata, const double timelim, int best_tour[tspdata->n]);
void two_opt(TSPdata* tspdata, double timelim, int best_tour[tspdata->n]);
void three_opt(TSPdata* tspdata, double timelim, int best_tour[tspdata->n]);

int compute_distance(double x1, double y1, double x2, double y2);
int compute_cost(TSPdata* tspdata, int* tour);
int is_feasible(TSPdata* tspdata, int* tour);

/***** open the file with given mode *****************************************/
FILE* open_file(char* fname, char* mode) {
  FILE* fp;
  fp = fopen(fname, mode);
  if (fp == NULL) {
    fprintf(stderr, "file not found: %s\n", fname);
    exit(EXIT_FAILURE);
  }
  return fp;
}

/***** malloc with error check ***********************************************/
void* malloc_e(size_t size) {
  void* s;
  if ((s = malloc(size)) == NULL) {
    fprintf(stderr, "malloc : not enough memory.\n");
    exit(EXIT_FAILURE);
  }
  return s;
}

/***** copy and read the parameters ******************************************/
/***** Feel free to modify this subroutine. **********************************/
void copy_parameters(int argc, char* argv[], Param* param) {
  /**** copy the default parameters ****/
  param->timelim = TIMELIM;
  param->givesol = GIVESOL;
  param->outformat = OUTFORMAT;
  strcpy(param->tourfile, TOURFILE);

  /**** read the parameters ****/
  if (argc > 0 && (argc % 2) == 0) {
    printf("USAGE: ./%s [param_name, param_value] [name, value]...\n", argv[0]);
    exit(EXIT_FAILURE);
  } else {
    int i;
    for (i = 1; i < argc; i += 2) {
      if (strcmp(argv[i], "timelim") == 0)
        param->timelim = atoi(argv[i + 1]);
      if (strcmp(argv[i], "givesol") == 0)
        param->givesol = atoi(argv[i + 1]);
      if (strcmp(argv[i], "outformat") == 0)
        param->outformat = atoi(argv[i + 1]);
      if (strcmp(argv[i], "tourfile") == 0)
        strcpy(param->tourfile, argv[i + 1]);
    }
  }
}

/***** prepare memory space **************************************************/
/***** Feel free to modify this subroutine. **********************************/
void prepare_memory(TSPdata* tspdata, Vdata* vdata) {
  int k, n;
  n = tspdata->n;
  tspdata->x = (double*)malloc_e(n * sizeof(double));
  tspdata->y = (double*)malloc_e(n * sizeof(double));
  vdata->bestsol = (int*)malloc_e(n * sizeof(int));
  /* the next line is just to give an initial solution */
  for (k = 0; k < n; k++)
    vdata->bestsol[k] = k;
}

/***** reading the header of a file in TSPLIB format *************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_header(FILE* in, TSPdata* tspdata) {
  char str[MAX_STR], name[MAX_STR], dim[MAX_STR], type[MAX_STR], edge[MAX_STR],
      min[MAX_STR];
  int flag = 0;

  for (;;) {
    char *w, *u;
    /* error */
    if (fgets(str, MAX_STR, in) == NULL) {
      fprintf(stderr, "error: invalid data input.\n");
      exit(EXIT_FAILURE);
    }
    /* halt condition */
    if (strcmp(str, "NODE_COORD_SECTION\n") == 0) {
      break;
    }
    if (strcmp(str, "TOUR_SECTION\n") == 0) {
      flag = 1;
      break;
    }
    /* data input */
    w = strtok(str, " :\n");
    u = strtok(NULL, " :\n");
    if (w == NULL || u == NULL)
      continue;
    if (strcmp("NAME", w) == 0)
      strcpy(name, u);
    if (strcmp("DIMENSION", w) == 0)
      strcpy(dim, u);
    if (strcmp("TYPE", w) == 0)
      strcpy(type, u);
    if (strcmp("EDGE_WEIGHT_TYPE", w) == 0)
      strcpy(edge, u);
    if (strcmp("MIN_NODE_NUM", w) == 0)
      strcpy(min, u);
  }

  /* read a TSP instance */
  if (flag == 0) {
    strcpy(tspdata->name, name);
    tspdata->min_node_num = atoi(min);
    tspdata->n = atoi(dim);
    if (strcmp("TSP", type) != 0 || strcmp("EUC_2D", edge) != 0) {
      fprintf(stderr, "error: invalid instance.\n");
      exit(EXIT_FAILURE);
    }
  }
  /* read a tour */
  else {
    if (strcmp("TOUR", type) != 0) {
      fprintf(stderr, "error: invalid tour.\n");
      exit(EXIT_FAILURE);
    }
  }
}

/***** reading the file of TSP instance **************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_tspfile(FILE* in, TSPdata* tspdata, Vdata* vdata) {
  char str[MAX_STR];
  int k;

  /* reading the instance */
  read_header(in, tspdata);
  prepare_memory(tspdata, vdata);
  for (k = 0; k < tspdata->n; k++) {
    int dummy;
    if (fgets(str, MAX_STR, in) == NULL)
      break;
    if (strcmp(str, "EOF\n") == 0)
      break;
    sscanf(str, "%d%lf%lf", &dummy, &(tspdata->x[k]), &(tspdata->y[k]));
  }
  if (k != tspdata->n) {
    fprintf(stderr, "error: invalid instance.\n");
    exit(EXIT_FAILURE);
  }
}

/***** read the tour in the TSPLIB format with feasibility check *************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_tourfile(FILE* in, TSPdata* tspdata, int* tour) {
  int k;

  read_header(in, tspdata);
  for (k = 0; k < tspdata->n; k++) {
    int val;
    if (fscanf(in, "%d", &val) == EOF)
      break;
    if (val == -1)
      break;
    tour[k] = val - 1;
  }
  if (k != tspdata->n) {
    fprintf(stderr, "error: invalid tour.\n");
    exit(EXIT_FAILURE);
  }
}

/***** output the tour in the TSPLIB format **********************************/
/***** note: the output tour starts from the node "1" ************************/
void output_tour(FILE* out, TSPdata* tspdata, int* tour) {
  int k, idx = 0;

  fprintf(out, "NAME : %s\n", tspdata->name);
  fprintf(out, "COMMENT : tour_length=%d\n", compute_cost(tspdata, tour));
  fprintf(out, "TYPE : TOUR\n");
  fprintf(out, "DIMENSION : %d\n", tspdata->n);
  fprintf(out, "TOUR_SECTION\n");
  for (k = 0; k < tspdata->n; k++)
    if (tour[k] == 0) {
      idx = k;
      break;
    }
  for (k = idx; k < tspdata->n; k++) {
    if (tour[k] < 0)
      break;
    fprintf(out, "%d\n", tour[k] + 1);
  }
  for (k = 0; k < idx; k++)
    fprintf(out, "%d\n", tour[k] + 1);
  fprintf(out, "-1\n");
  fprintf(out, "EOF\n");
}

/***** output the tour in the TSP_VIEW format
 * **********************************/
/***** note: the indices of the tour starts from "0"
 * ***************************/
void output_tour_for_tsp_view(FILE* out, TSPdata* tspdata, int* tour) {
  int k;

  fprintf(out, "%d\n", tspdata->n);
  for (k = 0; k < tspdata->n; k++) {
    fprintf(out, "%g %g\n", tspdata->x[k], tspdata->y[k]);
  }
  for (k = 0; k < tspdata->n; k++) {
    if (tour[k] < 0)
      break;
    fprintf(out, "%d\n", tour[k]);
  }
}

/***** check the feasibility and recompute the cost **************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void recompute_obj(Param* param, TSPdata* tspdata, Vdata* vdata) {
  if (!is_feasible(tspdata, vdata->bestsol)) {
    fprintf(stderr, "error: the computed tour is not feasible.\n");
    exit(EXIT_FAILURE);
  }
  printf("recomputed tour length = %d\n",
         compute_cost(tspdata, vdata->bestsol));
  printf("time for the search:   %7.2f seconds\n",
         vdata->endtime - vdata->starttime);
  printf("time to read the instance: %7.2f seconds\n",
         vdata->starttime - vdata->timebrid);
}

/***** cost of the tour ******************************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
int compute_cost(TSPdata* tspdata, int* tour) {
  int k, cost = 0, n;
  n = tspdata->n;

  for (k = 0; k < n - 1; k++) {
    if (tour[k + 1] < 0)
      break;
    cost += dist(tour[k], tour[k + 1]);
  }
  cost += dist(tour[k], tour[0]);
  return cost;
}

/***** check the feasibility of the tour *************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
int is_feasible(TSPdata* tspdata, int* tour) {
  int k, n, *visited, flag = 1, num_visited = 0;
  n = tspdata->n;
  visited = (int*)malloc_e(n * sizeof(int));

  for (k = 0; k < n; k++)
    visited[k] = 0;
  for (k = 0; k < n; k++) {
    if (tour[k] < 0) {
      break;
    }
    if (tour[k] >= n) {
      flag = 0;
      break;
    }
    if (visited[tour[k]]) {
      flag = 0;
      break;
    } else {
      visited[tour[k]] = 1;
      num_visited++;
    }
  }
  if (num_visited < tspdata->min_node_num)
    flag = 0;

  free(visited);
  /* if tour is feasible (resp., not feasible), then flag=1 (resp., 0) */
  return flag;
}

/***** sample algorithm ******************************************************/
void my_algorithm(Param* param, TSPdata* tspdata, Vdata* vdata) {
  for (int i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = -1;
  }

  nearest_neighbor(tspdata,
                   (param->timelim - cpu_time() + vdata->starttime) * 0.1,
                   vdata->bestsol);

  const double iter_tim_lim =
      (param->timelim - cpu_time() + vdata->starttime) / 20;

  while (cpu_time() - vdata->starttime < param->timelim) {
    if (tspdata->n > tspdata->min_node_num) {
      replace(tspdata,
              fmin(iter_tim_lim / 4,
                   param->timelim - cpu_time() + vdata->starttime),
              vdata->bestsol);
    }
    two_opt(
        tspdata,
        fmin(iter_tim_lim / 2, param->timelim - cpu_time() + vdata->starttime),
        vdata->bestsol);
    three_opt(
        tspdata,
        fmin(iter_tim_lim / 4, param->timelim - cpu_time() + vdata->starttime),
        vdata->bestsol);
  }
}

void nearest_neighbor(TSPdata* tspdata,
                      const double timelim,
                      int best_tour[tspdata->n]) {
  const int n_nodes = tspdata->n;
  const int n_min_nodes = tspdata->min_node_num;

  double starttime = cpu_time();
  int min_cost = INT_MAX;
  if (is_feasible(tspdata, best_tour)) {
    min_cost = compute_cost(tspdata, best_tour);
  }

  while (cpu_time() - starttime < timelim) {
    bool is_visiteds[n_nodes];
    int local_tour[n_nodes];
    for (int i = 0; i < n_nodes; i++) {
      is_visiteds[i] = false;
      local_tour[i] = -1;
    }

    const int first_node = rand() % n_nodes;
    local_tour[0] = first_node;
    is_visiteds[first_node] = true;

    for (int depth = 1; depth < n_min_nodes; depth++) {
      int nearest_cost = INT_MAX;

      for (int temp_next_node = 0; temp_next_node < n_nodes; temp_next_node++) {
        if (is_visiteds[temp_next_node]) {
          continue;
        }

        const int cost = dist(local_tour[depth - 1], temp_next_node);
        if (cost < nearest_cost) {
          nearest_cost = cost;
          local_tour[depth] = temp_next_node;
        }
      }
      is_visiteds[local_tour[depth]] = true;
    }

    const int cost = compute_cost(tspdata, local_tour);
    if (cost < min_cost) {
      min_cost = cost;
      for (int tour_idx = 0; tour_idx < n_nodes; tour_idx++) {
        best_tour[tour_idx] = local_tour[tour_idx];
      }
    }
  }
}

void replace(TSPdata* tspdata,
             const double timelim,
             int best_tour[tspdata->n]) {
  const int n_nodes = tspdata->n;
  const int n_min_nodes = tspdata->min_node_num;

  double starttime = cpu_time();
  int min_cost = INT_MAX;
  if (is_feasible(tspdata, best_tour)) {
    min_cost = compute_cost(tspdata, best_tour);
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
          dist(pre_insert, insert_node) + dist(insert_node, following_insert);
      const int delete_dist =
          dist(pre_delete, delete_node) + dist(delete_node, following_delete);

      reduced_dist = delete_dist - added_dist;
    } else {
      const int added_dist = dist(pre_insert, insert_node) +
                             dist(insert_node, following_insert) -
                             dist(pre_insert, following_insert);
      const int delete_dist = dist(pre_delete, delete_node) +
                              dist(delete_node, following_delete) -
                              dist(pre_delete, following_delete);

      reduced_dist = delete_dist - added_dist;
    }

    if (reduced_dist < 0) {
      continue;
    }

    int new_tour[n_min_nodes];
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
  }
}

void two_opt(TSPdata* tspdata, double timelim, int best_tour[tspdata->n]) {
  const int n_nodes = tspdata->n;
  const int n_min_nodes = tspdata->min_node_num;

  double starttime = cpu_time();

  while (cpu_time() - starttime < timelim) {
    const int a_tour_idx = rand() % n_min_nodes;
    const int c_tour_idx = (a_tour_idx + 1) % n_min_nodes;

    const int b_tour_idx =
        (a_tour_idx + 1 + rand() % (n_min_nodes - 1)) % n_min_nodes;
    const int d_tour_idx = (b_tour_idx + 1) % n_min_nodes;

    const int a_node = best_tour[a_tour_idx];
    const int b_node = best_tour[b_tour_idx];
    const int c_node = best_tour[c_tour_idx];
    const int d_node = best_tour[d_tour_idx];

    const int reduced_cost = dist(a_node, c_node) + dist(b_node, d_node) -
                             dist(a_node, b_node) - dist(c_node, d_node);
    if (reduced_cost <= 0) {
      continue;
    }

    int new_tour[n_min_nodes];
    int new_tour_idx = 0;
#if DEBUG
    for (int i = 0; i < n_min_nodes; i++) {
      new_tour[i] = -1;
    }
#endif

    for (int tour_idx = b_tour_idx; tour_idx != c_tour_idx;
         tour_idx = tour_idx == 0 ? n_min_nodes - 1 : (tour_idx - 1)) {
      new_tour[new_tour_idx] = best_tour[tour_idx];
      new_tour_idx++;
    }
    new_tour[new_tour_idx] = c_node;
    new_tour_idx++;

    for (int tour_idx = d_tour_idx; tour_idx != a_tour_idx;
         tour_idx = (tour_idx + 1) % n_min_nodes) {
      new_tour[new_tour_idx] = best_tour[tour_idx];
      new_tour_idx++;
    }
    new_tour[new_tour_idx] = a_node;
    // new_tour_idx++;

    for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
      best_tour[tour_idx] = new_tour[tour_idx];
    }
  }
}

void three_opt(TSPdata* tspdata, double timelim, int best_tour[tspdata->n]) {
  const int n_nodes = tspdata->n;
  const int n_min_nodes = tspdata->min_node_num;

  double starttime = cpu_time();

  while (cpu_time() - starttime < timelim) {
    int a_tour_idx, b_tour_idx, c_tour_idx;
    while (true) {
      a_tour_idx = rand() % n_min_nodes;
      b_tour_idx = rand() % n_min_nodes;
      c_tour_idx = rand() % n_min_nodes;

      if (a_tour_idx == b_tour_idx || b_tour_idx == c_tour_idx ||
          c_tour_idx == a_tour_idx) {
        continue;
      } else if ((a_tour_idx < b_tour_idx) + (b_tour_idx < c_tour_idx) +
                     (c_tour_idx < a_tour_idx) ==
                 1) {
        continue;
      }
      break;
    }
    const int d_tour_idx = (a_tour_idx + 1) % n_min_nodes;
    const int e_tour_idx = (b_tour_idx + 1) % n_min_nodes;
    const int f_tour_idx = (c_tour_idx + 1) % n_min_nodes;

    const int a_node = best_tour[a_tour_idx];
    const int d_node = best_tour[d_tour_idx];

    const int b_node = best_tour[b_tour_idx];
    const int e_node = best_tour[e_tour_idx];

    const int c_node = best_tour[c_tour_idx];
    const int f_node = best_tour[f_tour_idx];

    const int oririnal_cost =
        dist(a_node, d_node) + dist(b_node, e_node) + dist(c_node, f_node);

    int new_costs[4];
    new_costs[0] =
        dist(b_node, c_node) + dist(e_node, a_node) + dist(f_node, d_node);

    new_costs[1] =
        dist(b_node, a_node) + dist(f_node, e_node) + dist(c_node, d_node);

    new_costs[2] =
        dist(b_node, f_node) + dist(a_node, e_node) + dist(c_node, d_node);

    new_costs[3] =
        dist(b_node, f_node) + dist(a_node, c_node) + dist(e_node, d_node);

    int min_cost_idx;
    int min_new_cost = INT_MAX;
    for (int i = 0; i < 4; i++) {
      if (new_costs[i] < min_new_cost) {
        min_new_cost = new_costs[i];
        min_cost_idx = i;
      }
    }

    const int reduced_cost = oririnal_cost - min_new_cost;
    if (reduced_cost <= 0) {
      continue;
    }

    int new_tour[n_min_nodes];
    int new_tour_idx = 0;
#if DEBUG
    for (int i = 0; i < n_min_nodes; i++) {
      new_tour[i] = -1;
    }
#endif

    for (int tour_idx = d_tour_idx; tour_idx != b_tour_idx;
         tour_idx = (tour_idx + 1) % n_min_nodes) {
      new_tour[new_tour_idx] = best_tour[tour_idx];
      new_tour_idx++;
    }

    new_tour[new_tour_idx] = b_node;
    new_tour_idx++;

    if (min_cost_idx == 0) {
      for (int tour_idx = c_tour_idx; tour_idx != e_tour_idx;
           tour_idx = tour_idx == 0 ? n_min_nodes - 1 : (tour_idx - 1)) {
        new_tour[new_tour_idx] = best_tour[tour_idx];
        new_tour_idx++;
      }

      new_tour[new_tour_idx] = e_node;
      new_tour_idx++;

      for (int tour_idx = a_tour_idx; tour_idx != f_tour_idx;
           tour_idx = tour_idx == 0 ? n_min_nodes - 1 : (tour_idx - 1)) {
        new_tour[new_tour_idx] = best_tour[tour_idx];
        new_tour_idx++;
      }

      new_tour[new_tour_idx] = f_node;
      // new_tour_idx++;
    } else if (min_cost_idx == 1) {
      for (int tour_idx = a_tour_idx; tour_idx != f_tour_idx;
           tour_idx = tour_idx == 0 ? n_min_nodes - 1 : (tour_idx - 1)) {
        new_tour[new_tour_idx] = best_tour[tour_idx];
        new_tour_idx++;
      }
      new_tour[new_tour_idx] = f_node;
      new_tour_idx++;

      for (int tour_idx = e_tour_idx; tour_idx != c_tour_idx;
           tour_idx = (tour_idx + 1) % n_min_nodes) {
        new_tour[new_tour_idx] = best_tour[tour_idx];
        new_tour_idx++;
      }

      new_tour[new_tour_idx] = c_node;
      // new_tour_idx++;
    } else {
      for (int tour_idx = f_tour_idx; tour_idx != a_tour_idx;
           tour_idx = (tour_idx + 1) % n_min_nodes) {
        new_tour[new_tour_idx] = best_tour[tour_idx];
        new_tour_idx++;
      }

      new_tour[new_tour_idx] = a_node;
      new_tour_idx++;
      if (min_cost_idx == 2) {
        for (int tour_idx = e_tour_idx; tour_idx != c_tour_idx;
             tour_idx = (tour_idx + 1) % n_min_nodes) {
          new_tour[new_tour_idx] = best_tour[tour_idx];
          new_tour_idx++;
        }

        new_tour[new_tour_idx] = c_node;
        // new_tour_idx++;
      } else {
        for (int tour_idx = c_tour_idx; tour_idx != e_tour_idx;
             tour_idx = tour_idx == 0 ? n_min_nodes - 1 : (tour_idx - 1)) {
          new_tour[new_tour_idx] = best_tour[tour_idx];
          new_tour_idx++;
        }
        new_tour[new_tour_idx] = e_node;
        // new_tour_idx++;
      }
    }

    for (int tour_idx = 0; tour_idx < n_min_nodes; tour_idx++) {
      best_tour[tour_idx] = new_tour[tour_idx];
    }
  }
}

/***** main ******************************************************************/
int main(int argc, char* argv[]) {
  Param param;     /* parameters */
  TSPdata tspdata; /* data of TSP instance */
  Vdata vdata;     /* various data often needed during search */

  vdata.timebrid = cpu_time();
  copy_parameters(argc, argv, &param);
  read_tspfile(stdin, &tspdata, &vdata);
  if (param.givesol == 1)
    read_tourfile(stdin, &tspdata, vdata.bestsol);
  vdata.starttime = cpu_time();

  /*****

    Write your algorithm here.
    Of course you can add your subroutines outside main().
    At this point, the instance data is stored in the structure "tspdata".

       tspdata.name :     the name of the instance
       tspdata.n :        the number of nodes n
       tspdata.x[k] :     x-coordinate of node k (k = 0,1,...,n-1)
       tspdata.y[k] :     y-coordinate of node k (k = 0,1,...,n-1)
       tspdata.min_node_num :
                          the number of nodes that the solution must contain

   You can use the following macro to get the distance between node k and l.
       dist(k,l)    :  the distance between node k and l
                       (k,l = 0,1,...,n-1), where
                       dist(k,k) = 0 and
                       dist(k,l) = dist(l,k)
                       for any k and l

    Store your best tour (solution) in vdata.bestsol. If number of nodes in the
  tour is less than tspdata.n, you have to fill from bestsol[m] to bestsol[n-1]
  with -1. The "compute_cost()" will compute its objective value. The
  "is_feasible()" will return whether the solution is feasible or not. The
  format of vdata.bestsol is:

       vdata.bestsol[i] = k    if the i-th visited node of the tour is node k,
       vdata.bestsol[i] = -1   if i is greater than the number of visited nodes
  of the tour. where i,k = 0,1,...,n-1.

  *****/
  my_algorithm(&param, &tspdata, &vdata);

  vdata.endtime = cpu_time();
  recompute_obj(&param, &tspdata, &vdata);
  if (param.outformat == 1) {
    output_tour(open_file(param.tourfile, "w"), &tspdata, vdata.bestsol);
  } else if (param.outformat == 2) {
    output_tour_for_tsp_view(open_file(param.tourfile, "w"), &tspdata,
                             vdata.bestsol);
  }

  return EXIT_SUCCESS;
}
