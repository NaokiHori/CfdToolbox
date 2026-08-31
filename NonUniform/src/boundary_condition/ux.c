#include <stddef.h> // size_t
#include "boundary_condition.h"

// NOTE: implicit solver may need to be modified accordingly when the boundary conditions are changed

int impose_boundary_condition_ux_x(
    array_t * const ux_array
) {
  double (* const ux)[NX + 2] = ux_array->buffer;
  for (size_t j = 0; j <= NY + 1; j++) {
    ux[j][     0] = 0.;
    ux[j][     1] = 0.;
    ux[j][NX + 1] = 0.;
  }
  return 0;
}

int impose_boundary_condition_ux_y(
    array_t * const ux_array
) {
  const double ux_ym = 0.;
  const double ux_yp = 0.;
  double (* const ux)[NX + 2] = ux_array->buffer;
  for (size_t i = 0; i <= NX + 1; i++) {
    ux[     0][i] = ux_ym;
    ux[NY + 1][i] = ux_yp;
  }
  return 0;
}

