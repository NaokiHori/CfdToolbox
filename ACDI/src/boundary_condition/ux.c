#include <stddef.h> // size_t
#include "logger.h"
#include "boundary_condition.h"

int impose_boundary_condition_ux_x(
    const domain_t * const domain,
    double ** const ux
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t j = 0; j <= ny + 1; j++) {
    ux[j][     0] = 0.;
    ux[j][     1] = 0.;
    ux[j][nx + 1] = 0.;
  }
  return 0;
}

