#include <stddef.h>
#include "./impose_boundary_condition.h"

int impose_boundary_condition(
    const domain_t * const domain,
    double ** const concentration
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const hxxf = domain->hxxf;
  const double flux = 1.;
  for (size_t j = 1; j <= ny; j++) {
    concentration[     0][j] = concentration[1][j] + flux * hxxf[1];
    concentration[nx + 1][j] = 0.;
  }
  return 0;
}

