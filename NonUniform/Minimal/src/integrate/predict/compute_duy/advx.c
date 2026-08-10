#include "./advx.h"

int uy_advx(
    const domain_t * const domain,
    double ** const ux,
    double ** const uy,
    const double dt,
    double ** const duy
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  const double * const x_center_scale_factors = domain->x_coordinate.center_scale_factors;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
#pragma omp parallel for
  for (size_t j = 2; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double hx_c = x_center_scale_factors[i    ];
      const double hx_m = x_face_scale_factors[i    ];
      const double hx_p = x_face_scale_factors[i + 1];
      const double hy_c = y_face_scale_factors[j    ];
      const double hy_m = y_center_scale_factors[j - 1];
      const double hy_p = y_center_scale_factors[j    ];
      const double j_cc = hx_c * hy_c;
      const double j_mm = hx_m * hy_m;
      const double j_pm = hx_p * hy_m;
      const double j_mp = hx_m * hy_p;
      const double j_pp = hx_p * hy_p;
      const double ux_xm = + 0.5 * j_mm / hx_m * ux[j - 1][i    ]
                           + 0.5 * j_mp / hx_m * ux[j    ][i    ];
      const double ux_xp = + 0.5 * j_pm / hx_p * ux[j - 1][i + 1]
                           + 0.5 * j_pp / hx_p * ux[j    ][i + 1];
      const double duy_xm = - uy[j    ][i - 1]
                            + uy[j    ][i    ];
      const double duy_xp = - uy[j    ][i    ]
                            + uy[j    ][i + 1];
      duy[j][i] -= dt / j_cc * (
          + 0.5 * ux_xm * duy_xm
          + 0.5 * ux_xp * duy_xp
      );
    }
  }
  return 0;
}

