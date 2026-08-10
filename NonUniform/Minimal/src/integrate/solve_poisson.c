#include "logger.h"
#include "tridiagonal_solver.h"
#include "./solve_poisson.h"
#include "./transpose.h"

static int compute_matrix_vector_product(
    const size_t nitems,
    const double * const matrix,
    const double * const vector,
    double * const result
) {
  for (size_t j = 0; j < nitems; j++) {
    double * const value = result + j;
    *value = 0.;
    for (size_t i = 0; i < nitems; i++) {
      *value += matrix[j * nitems + i] * vector[i];
    }
  }
  return 0;
}

int solve_poisson(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const x_center_scale_factors = domain->x_coordinate.center_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  const double * const forward_operator = poisson_solver->forward_operator;
  const double * const backward_operator = poisson_solver->backward_operator;
  double ** const psi = flow_solver->psi;
  double * const buf0 = poisson_solver->buf0;
  double * const buf1 = poisson_solver->buf1;
  // assign right-hand side of Poisson equation
  {
    double ** const ux = flow_field->ux;
    double ** const uy = flow_field->uy;
    const double factor = 1. / dt;
#pragma omp parallel for
    for (size_t j = 1; j <= ny; j++) {
      for (size_t i = 1; i <= nx; i++) {
        const double dux = - ux[j    ][i    ]
                           + ux[j    ][i + 1];
        const double duy = - uy[j    ][i    ]
                           + uy[j + 1][i    ];
        const double div = (
            + 1. / x_center_scale_factors[i] * dux
            + 1. / y_center_scale_factors[j] * duy
        );
        buf0[(j - 1) * nx + (i - 1)] = factor * div;
      }
    }
  }
  // project x to spectral space
#pragma omp parallel for
  for (size_t j = 0; j < ny; j++) {
    compute_matrix_vector_product(nx, forward_operator, buf0 + j * nx, buf1 + j * nx);
  }
  // x-align to y-align
  if (0 != transpose(nx, ny, buf1, buf0)) {
    LOGGER_FAILURE("failed to transpose array from x-aligned to y-aligned");
    goto abort;
  }
  // solve linear systems in y
  {
    const double * const tridiagonal_solver_lower_diagonals = poisson_solver->tridiagonal_solver_lower_diagonals;
    const double * const tridiagonal_solver_main_diagonals = poisson_solver->tridiagonal_solver_main_diagonals;
    const double * const tridiagonal_solver_upper_diagonals = poisson_solver->tridiagonal_solver_upper_diagonals;
    double * const tridiagonal_solver_buffers = poisson_solver->tridiagonal_solver_buffers;
    const double * const eigenvalues = poisson_solver->eigenvalues;
#pragma omp parallel for
    for (size_t i = 0; i < nx; i++) {
      tridiagonal_solver_solve(
          ny,
          tridiagonal_solver_lower_diagonals,
          tridiagonal_solver_main_diagonals,
          tridiagonal_solver_upper_diagonals,
          eigenvalues[i],
          tridiagonal_solver_buffers + i * ny,
          buf0 + i * ny
      );
    }
  }
  // y-align to x-align
  if (0 != transpose(ny, nx, buf0, buf1)) {
    LOGGER_FAILURE("failed to transpose array from y-aligned to x-aligned");
    goto abort;
  }
  // project x to physical space
#pragma omp parallel for
  for (size_t j = 0; j < ny; j++) {
    compute_matrix_vector_product(nx, backward_operator, buf1 + j * nx, buf0 + j * nx);
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      psi[j][i] = buf0[(j - 1) * nx + (i - 1)];
    }
  }
  return 0;
abort:
  return 1;
}

