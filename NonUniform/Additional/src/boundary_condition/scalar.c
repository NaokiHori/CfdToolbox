#include <stddef.h> // size_t
#include "boundary_condition.h"

// NOTE: implicit solver may need to be modified accordingly when the boundary conditions are changed

int impose_boundary_condition_scalar_x(
    array_t * const scalar_array
) {
  double (* const scalar)[NX + 2] = scalar_array->buffer;
  for (size_t j = 0; j <= NY + 1; j++) {
    scalar[j][     0] = scalar[j][ 1];
    scalar[j][NX + 1] = scalar[j][NX];
  }
  return 0;
}

int impose_boundary_condition_scalar_y(
    array_t * const scalar_array
) {
  const double scalar_ym = - 0.5;
  const double scalar_yp = + 0.5;
  double (* const scalar)[NX + 2] = scalar_array->buffer;
  for (size_t i = 0; i <= NX + 1; i++) {
    scalar[     0][i] = scalar_ym;
    scalar[NY + 1][i] = scalar_yp;
  }
  return 0;
}

