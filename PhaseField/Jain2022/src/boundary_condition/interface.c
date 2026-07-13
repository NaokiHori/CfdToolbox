#include <stddef.h> // size_t
#include "logger.h"
#include "boundary_condition.h"

int impose_boundary_condition_vof_x(
    const domain_t * const domain,
    double ** const vof
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t j = 0; j <= ny + 1; j++) {
    vof[j][     0] = 0.;
    vof[j][nx + 1] = 0.;
  }
  return 0;
}

int impose_boundary_condition_normal_x_x(
    const domain_t * const domain,
    double ** const normal_x
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t j = 0; j <= ny + 1; j++) {
    normal_x[j][     0] = -1.;
    normal_x[j][nx + 1] = 1.;
  }
  return 0;
}

