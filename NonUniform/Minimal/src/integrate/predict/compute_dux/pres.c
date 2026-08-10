#include "./pres.h"

int ux_pres(
    const domain_t * const domain,
    double ** const p,
    const double dt,
    double ** const dux
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      dux[j][i] -= dt / x_face_scale_factors[i] * (
          - p[j    ][i - 1]
          + p[j    ][i    ]
      );
    }
  }
  return 0;
}

