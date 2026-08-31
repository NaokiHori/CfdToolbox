#include "./compute_dscalar.h"

int compute_dscalar(
    const double diffusivity,
    const domain_t * const domain,
    const array_t * const ux_array,
    const array_t * const uy_array,
    const array_t * const scalar_array,
    array_t * const dscalar_array,
    array_t * const dscalar_explicit_previous_array,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients
) {
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const coordinate_t * const y_coordinate = &domain->y_coordinate;
  const double * const x_face_scale_factors = x_coordinate->face_scale_factors;
  const double * const x_center_scale_factors = x_coordinate->center_scale_factors;
  const double * const y_face_scale_factors = y_coordinate->face_scale_factors;
  const double * const y_center_scale_factors = y_coordinate->center_scale_factors;
  const tridiagonal_t * const x_center_laplace_operators = x_coordinate->center_laplace_operators;
  const tridiagonal_t * const y_center_laplace_operators = y_coordinate->center_laplace_operators;
  const double (* restrict const ux)[NX + 2] = ux_array->buffer;
  const double (* restrict const uy)[NX + 2] = uy_array->buffer;
  const double (* restrict const scalar)[NX + 2] = scalar_array->buffer;
  double (* restrict const dscalar)[NX + 2] = dscalar_array->buffer;
  double (* restrict const dscalar_explicit_previous)[NX + 2] = dscalar_explicit_previous_array->buffer;
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    const double hy_c = y_center_scale_factors[j    ];
    const double hy_m = y_face_scale_factors[j    ];
    const double hy_p = y_face_scale_factors[j + 1];
    const tridiagonal_t * const y_center_laplace_operator = y_center_laplace_operators + j;
    const double y_lower_diagonal = y_center_laplace_operator->lower;
    const double y_main_diagonal = y_center_laplace_operator->main;
    const double y_upper_diagonal = y_center_laplace_operator->upper;
    for (size_t i = 1; i <= NX; i++) {
      const double hx_c = x_center_scale_factors[i    ];
      const double hx_m = x_face_scale_factors[i    ];
      const double hx_p = x_face_scale_factors[i + 1];
      const tridiagonal_t * const x_center_laplace_operator = x_center_laplace_operators + i;
      const double x_lower_diagonal = x_center_laplace_operator->lower;
      const double x_main_diagonal = x_center_laplace_operator->main;
      const double x_upper_diagonal = x_center_laplace_operator->upper;
      const double j_mc = hx_m * hy_c;
      const double j_cc = hx_c * hy_c;
      const double j_pc = hx_p * hy_c;
      const double j_cm = hx_c * hy_m;
      const double j_cp = hx_c * hy_p;
      const double ux_xm = j_mc / hx_m * ux[j    ][i    ];
      const double ux_xp = j_pc / hx_p * ux[j    ][i + 1];
      const double uy_ym = j_cm / hy_m * uy[j    ][i    ];
      const double uy_yp = j_cp / hy_p * uy[j + 1][i    ];
      const double dscalar_xm = - scalar[j    ][i - 1] + scalar[j    ][i    ];
      const double dscalar_xp = - scalar[j    ][i    ] + scalar[j    ][i + 1];
      const double dscalar_ym = - scalar[j - 1][i    ] + scalar[j    ][i    ];
      const double dscalar_yp = - scalar[j    ][i    ] + scalar[j + 1][i    ];
      const double advection_x = - 1. / j_cc * (
          + 0.5 * ux_xm * dscalar_xm
          + 0.5 * ux_xp * dscalar_xp
      );
      const double advection_y = - 1. / j_cc * (
          + 0.5 * uy_ym * dscalar_ym
          + 0.5 * uy_yp * dscalar_yp
      );
      const double diffusion_x = diffusivity * (
          + x_lower_diagonal * scalar[j    ][i - 1]
          + x_main_diagonal * scalar[j    ][i    ]
          + x_upper_diagonal * scalar[j    ][i + 1]
      );
      const double diffusion_y = diffusivity * (
          + y_lower_diagonal * scalar[j - 1][i    ]
          + y_main_diagonal * scalar[j    ][i    ]
          + y_upper_diagonal * scalar[j + 1][i    ]
      );
      dscalar[j][i] = dt * (
          + runge_kutta_coefficients->explicit_current * (advection_x + advection_y)
          + runge_kutta_coefficients->explicit_previous * dscalar_explicit_previous[j][i]
          + runge_kutta_coefficients->implicit * (diffusion_x + diffusion_y)
      );
      dscalar_explicit_previous[j][i] = advection_x + advection_y;
    }
  }
  return 0;
}

