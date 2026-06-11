#include <string.h>
#include "./halo.h"

static int copy(
    const double * const source,
    double * const destination
) {
  memcpy(destination, source, sizeof(double) * NDIRS);
  return 0;
}

int exchange_halo(
    df_t * const df
) {
  // periodic x boundaries
  for (size_t j = 0; j <= NY + 1; j++) {
    copy((*df)[j][NX], (*df)[j][     0]);
    copy((*df)[j][ 1], (*df)[j][NX + 1]);
  }
  return 0;
}

