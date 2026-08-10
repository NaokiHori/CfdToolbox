#include "./difx.h"

int uy_difx(
    const domain_t * const domain,
    const double c,
    double ** const uy,
    const double dt,
    double ** const duy
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  const double * const x_center_scale_factors = domain->x_coordinate.center_scale_factors;
#pragma omp parallel for
  for (size_t j = 2; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double lower_coefficient = 1. / x_center_scale_factors[i] / x_face_scale_factors[i    ];
      const double upper_coefficient = 1. / x_center_scale_factors[i] / x_face_scale_factors[i + 1];
      const double main_coefficient = - lower_coefficient - upper_coefficient;
      duy[j][i] += dt * c * (
          + lower_coefficient * uy[j    ][i - 1]
          + main_coefficient  * uy[j    ][i    ]
          + upper_coefficient * uy[j    ][i + 1]
      );
    }
  }
  return 0;
}

