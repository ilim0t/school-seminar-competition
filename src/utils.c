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

#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

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
void recompute_obj(const Param* const param,
                   const TSPdata* const tspdata,
                   const Vdata* const vdata) {
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
int compute_cost(const TSPdata* const tspdata, const int* const tour) {
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
int is_feasible(const TSPdata* const tspdata, const int* const tour) {
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
