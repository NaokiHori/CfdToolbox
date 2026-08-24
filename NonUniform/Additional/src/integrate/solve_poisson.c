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
  const double * const x_center_scale_factors = domain->x_coordinate.center_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  const double * const forward_operator = poisson_solver->forward_operator;
  const double * const backward_operator = poisson_solver->backward_operator;
  double (* const psi)[NX + 2] = flow_solver->psi->buffer;
  double * const buf0 = poisson_solver->buf0;
  double * const buf1 = poisson_solver->buf1;
  // assign right-hand side of Poisson equation
  {
    const array_t * const ux_array = flow_field->ux;
    const array_t * const uy_array = flow_field->uy;
    const double (* const ux)[NX + 2] = ux_array->buffer;
    const double (* const uy)[NX + 2] = uy_array->buffer;
    const double factor = 1. / dt;
#pragma omp parallel for
    for (size_t j = 1; j <= NY; j++) {
      for (size_t i = 1; i <= NX; i++) {
        const double dux = - ux[j    ][i    ]
                           + ux[j    ][i + 1];
        const double duy = - uy[j    ][i    ]
                           + uy[j + 1][i    ];
        const double div = (
            + 1. / x_center_scale_factors[i] * dux
            + 1. / y_center_scale_factors[j] * duy
        );
        buf0[(j - 1) * NX + (i - 1)] = factor * div;
      }
    }
  }
  // project x to spectral space
#pragma omp parallel for
  for (size_t j = 0; j < NY; j++) {
    compute_matrix_vector_product(NX, forward_operator, buf0 + j * NX, buf1 + j * NX);
  }
  // x-align to y-align
  if (0 != transpose(NX, NY, buf1, buf0)) {
    goto abort;
  }
  // solve linear systems in y
  {
    const tridiagonal_t * const diagonals = poisson_solver->diagonals;
    double * const tridiagonal_solver_buffers = poisson_solver->tridiagonal_solver_buffers;
    const double * const eigenvalues = poisson_solver->eigenvalues;
#pragma omp parallel for
    for (size_t i = 0; i < NX; i++) {
      tridiagonal_solver_solve(
          NY,
          diagonals,
          1.,
          eigenvalues[i],
          tridiagonal_solver_buffers + i * NY,
          buf0 + i * NY
      );
    }
  }
  // y-align to x-align
  if (0 != transpose(NY, NX, buf0, buf1)) {
    goto abort;
  }
  // project x to physical space
#pragma omp parallel for
  for (size_t j = 0; j < NY; j++) {
    compute_matrix_vector_product(NX, backward_operator, buf1 + j * NX, buf0 + j * NX);
  }
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      psi[j][i] = buf0[(j - 1) * NX + (i - 1)];
    }
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

