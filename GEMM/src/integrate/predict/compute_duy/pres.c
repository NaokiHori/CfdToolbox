#include "./pres.h"

int uy_pres(
    const domain_t * const domain,
    double ** const p,
    const double dt,
    double ** const duy
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
#pragma omp parallel for
  for (size_t j = 2; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      duy[j][i] -= dt / y_face_scale_factors[j] * (
          - p[j - 1][i    ]
          + p[j    ][i    ]
      );
    }
  }
  return 0;
}

