#include <stddef.h>
#include "./exchange_halo.h"

int exchange_halo(
    const domain_t * const domain,
    double ** const concentration
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t i = 0; i <= nx + 1; i++) {
    concentration[i][     0] = concentration[i][ny];
    concentration[i][ny + 1] = concentration[i][ 1];
  }
  return 0;
}
