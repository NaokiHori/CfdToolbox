#include <stddef.h> // size_t
#include "boundary_condition.h"

// NOTE: implicit solver may need to be modified accordingly when the boundary conditions are changed

int impose_boundary_condition_uy_x(
    array_t * const uy_array
) {
  const double uy_xm = 0.;
  const double uy_xp = 0.;
  double (* const uy)[NX + 2] = uy_array->buffer;
  for (size_t j = 0; j <= NY + 1; j++) {
    uy[j][     0] = uy_xm;
    uy[j][NX + 1] = uy_xp;
  }
  return 0;
}

int impose_boundary_condition_uy_y(
    array_t * const uy_array
) {
  double (* const uy)[NX + 2] = uy_array->buffer;
  for (size_t i = 0; i <= NX + 1; i++) {
    uy[     0][i] = 0.;
    uy[     1][i] = 0.;
    uy[NY + 1][i] = 0.;
  }
  return 0;
}

