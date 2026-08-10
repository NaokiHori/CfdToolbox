#include "./dify.h"

int ux_dify(
    const domain_t * const domain,
    const double c,
    double ** const ux,
    const double dt,
    double ** const dux
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    const double lower_coefficient = 1. / y_center_scale_factors[j] / y_face_scale_factors[j    ];
    const double upper_coefficient = 1. / y_center_scale_factors[j] / y_face_scale_factors[j + 1];
    const double main_coefficient = - lower_coefficient - upper_coefficient;
    for (size_t i = 2; i <= nx; i++) {
      dux[j][i] += dt * c * (
          + lower_coefficient * ux[j - 1][i    ]
          + main_coefficient  * ux[j    ][i    ]
          + upper_coefficient * ux[j + 1][i    ]
      );
    }
  }
  return 0;
}

