#include <stddef.h> // size_t
#include "boundary_condition.h"

int impose_boundary_condition_uy_x(
    const domain_t * const domain,
    double ** const uy
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double uy_xm = 0.;
  const double uy_xp = 0.;
  for (size_t j = 0; j <= ny + 1; j++) {
    uy[j][     0] = uy_xm;
    uy[j][nx + 1] = uy_xp;
  }
  return 0;
}

int impose_boundary_condition_uy_y(
    const domain_t * const domain,
    double ** const uy
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t i = 0; i <= nx + 1; i++) {
    uy[     0][i] = 0.;
    uy[     1][i] = 0.;
    uy[ny + 1][i] = 0.;
  }
  return 0;
}

