#include "./compute_dux.h"

int compute_dux(
    const parameters_t * const parameters,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const flow_solver_t * const flow_solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients
) {
  const double diffusivity = parameters->fluid_diffusivity;
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const coordinate_t * const y_coordinate = &domain->y_coordinate;
  const double * const x_face_scale_factors = x_coordinate->face_scale_factors;
  const double * const x_center_scale_factors = x_coordinate->center_scale_factors;
  const double * const y_face_scale_factors = y_coordinate->face_scale_factors;
  const double * const y_center_scale_factors = y_coordinate->center_scale_factors;
  const tridiagonal_t * const x_face_laplace_operators = x_coordinate->face_laplace_operators;
  const tridiagonal_t * const y_center_laplace_operators = y_coordinate->center_laplace_operators;
  const array_t * const ux_array = flow_field->ux;
  const array_t * const uy_array = flow_field->uy;
  const array_t * const p_array = flow_field->p;
  const double (* restrict const ux)[NX + 2] = ux_array->buffer;
  const double (* restrict const uy)[NX + 2] = uy_array->buffer;
  const double (* restrict const p)[NX + 2] = p_array->buffer;
  double (* restrict const dux)[NX + 2] = flow_solver->dux->buffer;
  double (* restrict const dux_explicit_previous)[NX + 2] = flow_solver->dux_explicit_previous->buffer;
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    const double hy_c = y_center_scale_factors[j    ];
    const double hy_m = y_face_scale_factors[j    ];
    const double hy_p = y_face_scale_factors[j + 1];
    const tridiagonal_t * const y_center_laplace_operator = y_center_laplace_operators + j;
    const double y_lower_diagonal = y_center_laplace_operator->lower;
    const double y_main_diagonal = y_center_laplace_operator->main;
    const double y_upper_diagonal = y_center_laplace_operator->upper;
    for (size_t i = 2; i <= NX; i++) {
      const double hx_c = x_face_scale_factors[i    ];
      const double hx_m = x_center_scale_factors[i - 1];
      const double hx_p = x_center_scale_factors[i    ];
      const tridiagonal_t * const x_face_laplace_operator = x_face_laplace_operators + i;
      const double x_lower_diagonal = x_face_laplace_operator->lower;
      const double x_main_diagonal = x_face_laplace_operator->main;
      const double x_upper_diagonal = x_face_laplace_operator->upper;
      const double j_mm = hx_m * hy_m;
      const double j_mc = hx_m * hy_c;
      const double j_mp = hx_m * hy_p;
      const double j_cc = hx_c * hy_c;
      const double j_pm = hx_p * hy_m;
      const double j_pc = hx_p * hy_c;
      const double j_pp = hx_p * hy_p;
      const double ux_xm = + 0.5 * j_mc / hx_m * ux[j    ][i - 1]
                           + 0.5 * j_cc / hx_c * ux[j    ][i    ];
      const double ux_xp = + 0.5 * j_cc / hx_c * ux[j    ][i    ]
                           + 0.5 * j_pc / hx_p * ux[j    ][i + 1];
      const double uy_ym = + 0.5 * j_mm / hy_m * uy[j    ][i - 1]
                           + 0.5 * j_pm / hy_m * uy[j    ][i    ];
      const double uy_yp = + 0.5 * j_mp / hy_p * uy[j + 1][i - 1]
                           + 0.5 * j_pp / hy_p * uy[j + 1][i    ];
      const double dux_xm = - ux[j    ][i - 1] + ux[j    ][i    ];
      const double dux_xp = - ux[j    ][i    ] + ux[j    ][i + 1];
      const double dux_ym = - ux[j - 1][i    ] + ux[j    ][i    ];
      const double dux_yp = - ux[j    ][i    ] + ux[j + 1][i    ];
      const double advection_x = - 1. / j_cc * (
          + 0.5 * ux_xm * dux_xm
          + 0.5 * ux_xp * dux_xp
      );
      const double advection_y = - 1. / j_cc * (
          + 0.5 * uy_ym * dux_ym
          + 0.5 * uy_yp * dux_yp
      );
      const double diffusion_x = diffusivity * (
          + x_lower_diagonal * ux[j    ][i - 1]
          + x_main_diagonal * ux[j    ][i    ]
          + x_upper_diagonal * ux[j    ][i + 1]
      );
      const double diffusion_y = diffusivity * (
          + y_lower_diagonal * ux[j - 1][i    ]
          + y_main_diagonal * ux[j    ][i    ]
          + y_upper_diagonal * ux[j + 1][i    ]
      );
      const double pressure = - 1. / hx_c * (
          - p[j    ][i - 1]
          + p[j    ][i    ]
      );
      dux[j][i] = dt * (
          + runge_kutta_coefficients->explicit_current * (advection_x + advection_y)
          + runge_kutta_coefficients->explicit_previous * dux_explicit_previous[j][i]
          + runge_kutta_coefficients->implicit * (diffusion_x + diffusion_y + pressure)
      );
      dux_explicit_previous[j][i] = advection_x + advection_y;
    }
  }
  return 0;
}

