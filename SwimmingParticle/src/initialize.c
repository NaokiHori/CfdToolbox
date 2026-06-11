#include <stdlib.h>
#include "constant.h"
#include "./exchange_halo.h"
#include "./impose_boundary_condition.h"
#include "./initialize.h"

int initialize(
    const domain_t * const domain,
    double ** const concentration
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t j = 1; j <= ny; j++) {
      concentration[i][j] = 1e-8 * rand() / RAND_MAX;
    }
  }
  impose_boundary_condition(domain, concentration);
  exchange_halo(domain, concentration);
  return 0;
}
