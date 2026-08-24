#include "./compute_duy.h"

int compute_duy(
    const parameters_t * const parameters,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const flow_solver_t * const flow_solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients
) {
  const double density_ratio = parameters->density_ratio;
  const double diffusivity = parameters->fluid_diffusivity;
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const coordinate_t * const y_coordinate = &domain->y_coordinate;
  const double * const x_face_scale_factors = x_coordinate->face_scale_factors;
  const double * const x_center_scale_factors = x_coordinate->center_scale_factors;
  const double * const y_face_scale_factors = y_coordinate->face_scale_factors;
  const double * const y_center_scale_factors = y_coordinate->center_scale_factors;
  const tridiagonal_t * const x_center_laplace_operators = x_coordinate->center_laplace_operators;
  const tridiagonal_t * const y_face_laplace_operators = y_coordinate->face_laplace_operators;
  const array_t * const ux_array = flow_field->ux;
  const array_t * const uy_array = flow_field->uy;
  const array_t * const p_array = flow_field->p;
  const array_t * const temperature_array = flow_field->temperature;
  const array_t * const salinity_array = flow_field->salinity;
  const double (* restrict const ux)[NX + 2] = ux_array->buffer;
  const double (* restrict const uy)[NX + 2] = uy_array->buffer;
  const double (* restrict const p)[NX + 2] = p_array->buffer;
  const double (* restrict const temperature)[NX + 2] = temperature_array->buffer;
  const double (* restrict const salinity)[NX + 2] = salinity_array->buffer;
  double (* restrict const duy)[NX + 2] = flow_solver->duy->buffer;
  double (* restrict const duy_explicit_previous)[NX + 2] = flow_solver->duy_explicit_previous->buffer;
#pragma omp parallel for
  for (size_t j = 2; j <= NY; j++) {
    const double hy_c = y_face_scale_factors[j    ];
    const double hy_m = y_center_scale_factors[j - 1];
    const double hy_p = y_center_scale_factors[j    ];
    const tridiagonal_t * const y_face_laplace_operator = y_face_laplace_operators + j;
    const double y_lower_diagonal = y_face_laplace_operator->lower;
    const double y_main_diagonal = y_face_laplace_operator->main;
    const double y_upper_diagonal = y_face_laplace_operator->upper;
    for (size_t i = 1; i <= NX; i++) {
      const double hx_c = x_center_scale_factors[i    ];
      const double hx_m = x_face_scale_factors[i    ];
      const double hx_p = x_face_scale_factors[i + 1];
      const tridiagonal_t * const x_center_laplace_operator = x_center_laplace_operators + i;
      const double x_lower_diagonal = x_center_laplace_operator->lower;
      const double x_main_diagonal = x_center_laplace_operator->main;
      const double x_upper_diagonal = x_center_laplace_operator->upper;
      const double j_mm = hx_m * hy_m;
      const double j_mp = hx_m * hy_p;
      const double j_cm = hx_c * hy_m;
      const double j_cc = hx_c * hy_c;
      const double j_cp = hx_c * hy_p;
      const double j_pm = hx_p * hy_m;
      const double j_pp = hx_p * hy_p;
      const double ux_xm = + 0.5 * j_mm / hx_m * ux[j - 1][i    ]
                           + 0.5 * j_mp / hx_m * ux[j    ][i    ];
      const double ux_xp = + 0.5 * j_pm / hx_p * ux[j - 1][i + 1]
                           + 0.5 * j_pp / hx_p * ux[j    ][i + 1];
      const double uy_ym = + 0.5 * j_cm / hy_m * uy[j - 1][i    ]
                           + 0.5 * j_cc / hy_c * uy[j    ][i    ];
      const double uy_yp = + 0.5 * j_cc / hy_c * uy[j    ][i    ]
                           + 0.5 * j_cp / hy_p * uy[j + 1][i    ];
      const double duy_xm = - uy[j    ][i - 1] + uy[j    ][i    ];
      const double duy_xp = - uy[j    ][i    ] + uy[j    ][i + 1];
      const double duy_ym = - uy[j - 1][i    ] + uy[j    ][i    ];
      const double duy_yp = - uy[j    ][i    ] + uy[j + 1][i    ];
      const double advection_x = - 1. / j_cc * (
          + 0.5 * ux_xm * duy_xm
          + 0.5 * ux_xp * duy_xp
      );
      const double advection_y = - 1. / j_cc * (
          + 0.5 * uy_ym * duy_ym
          + 0.5 * uy_yp * duy_yp
      );
      const double pressure = - 1. / hy_c * (
          - p[j - 1][i    ]
          + p[j    ][i    ]
      );
      const double diffusion_x = diffusivity * (
          + x_lower_diagonal * uy[j    ][i - 1]
          + x_main_diagonal * uy[j    ][i    ]
          + x_upper_diagonal * uy[j    ][i + 1]
      );
      const double diffusion_y = diffusivity * (
          + y_lower_diagonal * uy[j - 1][i    ]
          + y_main_diagonal * uy[j    ][i    ]
          + y_upper_diagonal * uy[j + 1][i    ]
      );
      const double buoyancy =
        (0.5 * temperature[j - 1][i] + 0.5 * temperature[j][i])
        -
        density_ratio * (0.5 * salinity[j - 1][i] + 0.5 * salinity[j][i]);
      duy[j][i] = dt * (
          + runge_kutta_coefficients->explicit_current * (advection_x + advection_y + buoyancy)
          + runge_kutta_coefficients->explicit_previous * duy_explicit_previous[j][i]
          + runge_kutta_coefficients->implicit * (diffusion_x + diffusion_y + pressure)
      );
      duy_explicit_previous[j][i] = advection_x + advection_y + buoyancy;
    }
  }
  return 0;
}

