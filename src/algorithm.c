#include "algorithm.h"
#include "cpu_time.h"

void my_algorithm(Param* param, TSPdata* tspdata, Vdata* vdata) {
  int i, tmp;

  for (i = 0; i < tspdata->n; i++) {
    vdata->bestsol[i] = i;
  }
  while (cpu_time() - vdata->starttime < param->timelim) {
    for (i = 0; i < tspdata->n; i++) {
      tmp = vdata->bestsol[i];
      vdata->bestsol[i] = vdata->bestsol[(i + 1) % (tspdata->n)];
      vdata->bestsol[(i + 1) % (tspdata->n)] = tmp;
    }
  }
}
