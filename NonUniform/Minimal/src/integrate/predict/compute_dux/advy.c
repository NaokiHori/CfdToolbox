#include "./advy.h"

int ux_advy(
    const domain_t * const domain,
    double ** const uy,
    double ** const ux,
    const double dt,
    double ** const dux
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  const double * const x_center_scale_factors = domain->x_coordinate.center_scale_factors;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      const double hx_c = x_face_scale_factors[i    ];
      const double hx_m = x_center_scale_factors[i - 1];
      const double hx_p = x_center_scale_factors[i    ];
      const double hy_c = y_center_scale_factors[j    ];
      const double hy_m = y_face_scale_factors[j    ];
      const double hy_p = y_face_scale_factors[j + 1];
      const double j_cc = hx_c * hy_c;
      const double j_mm = hx_m * hy_m;
      const double j_pm = hx_p * hy_m;
      const double j_mp = hx_m * hy_p;
      const double j_pp = hx_p * hy_p;
      const double uy_ym = + 0.5 * j_mm / hy_m * uy[j    ][i - 1]
                           + 0.5 * j_pm / hy_m * uy[j    ][i    ];
      const double uy_yp = + 0.5 * j_mp / hy_p * uy[j + 1][i - 1]
                           + 0.5 * j_pp / hy_p * uy[j + 1][i    ];
      const double dux_ym = - ux[j - 1][i    ]
                            + ux[j    ][i    ];
      const double dux_yp = - ux[j    ][i    ]
                            + ux[j + 1][i    ];
      dux[j][i] -= dt / j_cc * (
          + 0.5 * uy_ym * dux_ym
          + 0.5 * uy_yp * dux_yp
      );
    }
  }
  return 0;
}

