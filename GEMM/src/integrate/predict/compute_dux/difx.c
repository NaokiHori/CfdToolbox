#include "./difx.h"

int ux_difx(
    const domain_t * const domain,
    const double c,
    double ** const ux,
    const double dt,
    double ** const dux
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  const double * const x_center_scale_factors = domain->x_coordinate.center_scale_factors;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      const double lower_coefficient = 1. / x_face_scale_factors[i] / x_center_scale_factors[i - 1];
      const double upper_coefficient = 1. / x_face_scale_factors[i] / x_center_scale_factors[i    ];
      const double main_coefficient = - lower_coefficient - upper_coefficient;
      dux[j][i] += dt * c * (
          + lower_coefficient * ux[j    ][i - 1]
          + main_coefficient  * ux[j    ][i    ]
          + upper_coefficient * ux[j    ][i + 1]
      );
    }
  }
  return 0;
}

