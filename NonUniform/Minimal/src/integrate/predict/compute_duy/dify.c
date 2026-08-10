#include "./dify.h"

int uy_dify(
    const domain_t * const domain,
    const double c,
    double ** const uy,
    const double dt,
    double ** const duy
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
#pragma omp parallel for
  for (size_t j = 2; j <= ny; j++) {
    const double lower_coefficient = 1. / y_face_scale_factors[j] / y_center_scale_factors[j - 1];
    const double upper_coefficient = 1. / y_face_scale_factors[j] / y_center_scale_factors[j    ];
    const double main_coefficient = - lower_coefficient - upper_coefficient;
    for (size_t i = 1; i <= nx; i++) {
      duy[j][i] += dt * c * (
          + lower_coefficient * uy[j - 1][i    ]
          + main_coefficient  * uy[j    ][i    ]
          + upper_coefficient * uy[j + 1][i    ]
      );
    }
  }
  return 0;
}

