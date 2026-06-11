#include <stddef.h>
#include "./stream.h"

static int swap(
    double * const a,
    double * const b
) {
  double tmp = *a;
  *a = *b;
  *b = tmp;
  return 0;
}

int process_streaming_inplace(
    df_t * const df
) {
#pragma omp parallel for
  for (size_t j = 0; j <= NY; j++) {
    swap((*df)[j][0] + 1, (*df)[j    ][1] + 3);
    swap((*df)[j][0] + 2, (*df)[j + 1][0] + 4);
    swap((*df)[j][0] + 5, (*df)[j + 1][1] + 7);
    for (size_t i = 1; i <= NX; i++) {
      swap((*df)[j][i] + 1, (*df)[j    ][i + 1] + 3);
      swap((*df)[j][i] + 2, (*df)[j + 1][i    ] + 4);
      swap((*df)[j][i] + 5, (*df)[j + 1][i + 1] + 7);
      swap((*df)[j][i] + 6, (*df)[j + 1][i - 1] + 8);
    }
    swap((*df)[j][NX + 1] + 6, (*df)[j + 1][NX] + 8);
  }
  for (size_t i = 0; i <= NX; i++) {
    swap((*df)[NY + 1][i] + 1, (*df)[NY + 1][i + 1] + 3);
  }
#pragma omp parallel for
  for (size_t j = 0; j <= NY + 1; j++) {
    for (size_t i = 1; i <= NX; i++) {
      double * const dfij = (*df)[j][i];
      swap(dfij + 1, dfij + 3);
      swap(dfij + 2, dfij + 4);
      swap(dfij + 5, dfij + 7);
      swap(dfij + 6, dfij + 8);
    }
  }
  return 0;
}

