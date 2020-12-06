#include "algorithm.h"
#include "cpu_time.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

void swap(int a[3], int b[3]) {
  for (int i = 0; i < 3; i++) {
    int temp = a[i];
    a[i] = b[i];
    b[i] = temp;
  }
}

void shaker_sort(int n, int array[n][3]) {
  for (int i = 0; i < n; i++) {
    // head to tail
    for (int j = i + 1; j < n; j++) {
      if (array[j][0] < array[j - 1][0]) {
        swap(array[j], array[j - 1]);
      }
    }
    n--;

    // tail to head
    for (int k = n - 1; k > i; k--) {
      if (array[k][0] < array[k - 1][0]) {
        swap(array[k], array[k - 1]);
      }
    }
  }
}

void print_tour(const int* const local_tour, const TSPdata* const tspdata) {
  printf("%d:\t", compute_cost(tspdata, local_tour));
  for (int i = 0; i < tspdata->n; i++) {
    if (local_tour[i] < 0) {
      break;
    }
    printf("%d, ", local_tour[i]);
  }

  printf("\n");
}