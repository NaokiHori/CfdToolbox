#include "./advy.h"

int uy_advy(
    const domain_t * const domain,
    double ** const uy,
    const double dt,
    double ** const duy
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const x_center_scale_factors = domain->x_coordinate.center_scale_factors;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
#pragma omp parallel for
  for (size_t j = 2; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double hx = x_center_scale_factors[i];
      const double hy_c = y_face_scale_factors[j];
      const double hy_m = y_center_scale_factors[j - 1];
      const double hy_p = y_center_scale_factors[j    ];
      const double j_m = hx * hy_m;
      const double j_c = hx * hy_c;
      const double j_p = hx * hy_p;
      const double uy_ym = + 0.5 * j_m / hy_m * uy[j - 1][i    ]
                           + 0.5 * j_c / hy_c * uy[j    ][i    ];
      const double uy_yp = + 0.5 * j_c / hy_c * uy[j    ][i    ]
                           + 0.5 * j_p / hy_p * uy[j + 1][i    ];
      const double duy_ym = - uy[j - 1][i    ]
                            + uy[j    ][i    ];
      const double duy_yp = - uy[j    ][i    ]
                            + uy[j + 1][i    ];
      duy[j][i] -= dt / j_c * (
          + 0.5 * uy_ym * duy_ym
          + 0.5 * uy_yp * duy_yp
      );
    }
  }
  return 0;
}

