#include "algorithm.h"
#include "cpu_time.h"
#include "utils.h"

#include <stdlib.h>
#include <time.h>

int main(int argc, char* argv[]) {
  Param param;     /* parameters */
  TSPdata tspdata; /* data of TSP instance */
  Vdata vdata;     /* various data often needed during search */

  vdata.timebrid = cpu_time();
  copy_parameters(argc, argv, &param);
  FILE* fp = fopen("../data/test.tsp", "r");
  read_tspfile(fp, &tspdata, &vdata);
  if (param.givesol == 1)
    read_tourfile(fp, &tspdata, vdata.bestsol);
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
  struct timespec ts;
  if (timespec_get(&ts, TIME_UTC) == 0) {
    printf("error");
  } else {
    srandom(ts.tv_nsec ^ ts.tv_sec);
  }

  // nearest_neighbor_algorithm(&param, &tspdata, &vdata);
  // insertion_algorithm(&param, &tspdata, &vdata);
  // two_approximation_algorithm(&param, &tspdata, &vdata);

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
