#pragma once
#include <stdio.h>

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
#define TIMELIM 3 /* the time limit for the algorithm in seconds */
#define GIVESOL 0 /* 1: input a solution; 0: do not give a solution */
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

int compute_distance(double x1, double y1, double x2, double y2);
int compute_cost(TSPdata* tspdata, int* tour);
int is_feasible(TSPdata* tspdata, int* tour);