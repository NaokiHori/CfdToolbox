#include "./update_pressure.h"

int update_pressure(
    const parameters_t * const parameters,
    const domain_t * const domain,
    const double dt,
    const flow_solver_t * const flow_solver,
    flow_field_t * const flow_field
) {
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const coordinate_t * const y_coordinate = &domain->y_coordinate;
  const tridiagonal_t * const x_center_laplace_operators = x_coordinate->center_laplace_operators;
  const tridiagonal_t * const y_center_laplace_operators = y_coordinate->center_laplace_operators;
  const array_t * const psi_array = flow_solver->psi;
  array_t * const p_array = flow_field->p;
  const double (* const psi)[NX + 2] = psi_array->buffer;
  double (* const p)[NX + 2] = p_array->buffer;
  const double diffusivity = parameters->fluid_diffusivity;
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    const tridiagonal_t * const y_center_laplace_operator = y_center_laplace_operators + j;
    const double y_lower_diagonal = y_center_laplace_operator->lower;
    const double y_main_diagonal = y_center_laplace_operator->main;
    const double y_upper_diagonal = y_center_laplace_operator->upper;
    for (size_t i = 1; i <= NX; i++) {
      const tridiagonal_t * const x_center_laplace_operator = x_center_laplace_operators + i;
      const double x_lower_diagonal = x_center_laplace_operator->lower;
      const double x_main_diagonal = x_center_laplace_operator->main;
      const double x_upper_diagonal = x_center_laplace_operator->upper;
      p[j][i] +=
        psi[j][i]
        -
        0.5 * dt * diffusivity * (
            + x_lower_diagonal * psi[j    ][i - 1]
            + x_main_diagonal * psi[j    ][i    ]
            + x_upper_diagonal * psi[j    ][i + 1]
            + y_lower_diagonal * psi[j - 1][i    ]
            + y_main_diagonal * psi[j    ][i    ]
            + y_upper_diagonal * psi[j + 1][i    ]
        );
    }
  }
  return 0;
}

