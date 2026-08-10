#include <stddef.h> // size_t
#include "boundary_condition.h"

int impose_boundary_condition_ux_x(
    const domain_t * const domain,
    double ** const ux
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  for (size_t j = 0; j <= ny + 1; j++) {
    ux[j][     0] = 0.;
    ux[j][     1] = 0.;
    ux[j][nx + 1] = 0.;
  }
  return 0;
}

int impose_boundary_condition_ux_y(
    const domain_t * const domain,
    double ** const ux
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double ux_ym = 0.;
  const double ux_yp = 1.;
  for (size_t i = 0; i <= nx + 1; i++) {
    ux[     0][i] = ux_ym;
    ux[ny + 1][i] = ux_yp;
  }
  return 0;
}

