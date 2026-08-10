#include "./advx.h"

int ux_advx(
    const domain_t * const domain,
    double ** const ux,
    const double dt,
    double ** const dux
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  const double * const x_center_scale_factors = domain->x_coordinate.center_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      const double hx_c = x_face_scale_factors[i];
      const double hx_m = x_center_scale_factors[i - 1];
      const double hx_p = x_center_scale_factors[i    ];
      const double hy = y_center_scale_factors[j];
      const double j_m = hx_m * hy;
      const double j_c = hx_c * hy;
      const double j_p = hx_p * hy;
      const double ux_xm = + 0.5 * j_m / hx_m * ux[j    ][i - 1]
                           + 0.5 * j_c / hx_c * ux[j    ][i    ];
      const double ux_xp = + 0.5 * j_c / hx_c * ux[j    ][i    ]
                           + 0.5 * j_p / hx_p * ux[j    ][i + 1];
      const double dux_xm = - ux[j    ][i - 1]
                            + ux[j    ][i    ];
      const double dux_xp = - ux[j    ][i    ]
                            + ux[j    ][i + 1];
      dux[j][i] -= dt / j_c * (
          + 0.5 * ux_xm * dux_xm
          + 0.5 * ux_xp * dux_xp
      );
    }
  }
  return 0;
}

