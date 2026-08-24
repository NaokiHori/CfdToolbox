#include "boundary_condition.h"
#include "logger.h"
#include "tridiagonal_solver.h"
#include "./update_scalar.h"

int update_scalar(
    const double diffusivity,
    const helmholtz_solver_t * const solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients,
    array_t * const dscalar_array,
    array_t * const scalar_array
) {
  double * const buffers = solver->buffers;
  double * const transposed = solver->transposed;
  const tridiagonal_t * const x_diagonals = solver->x_diagonals;
  const tridiagonal_t * const y_diagonals = solver->y_diagonals;
  double (* const dscalar)[NX + 2] = dscalar_array->buffer;
  double (* const scalar)[NX + 2] = scalar_array->buffer;
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    tridiagonal_solver_solve(
        NX,
        x_diagonals,
        - 0.5 * diffusivity * runge_kutta_coefficients->implicit * dt,
        1.,
        buffers + (j - 1) * NX,
        &dscalar[j][1]
    );
  }
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      transposed[(i - 1) * NY + (j - 1)] = dscalar[j][i];
    }
  }
#pragma omp parallel for
  for (size_t i = 0; i < NX; i++) {
    tridiagonal_solver_solve(
        NY,
        y_diagonals,
        - 0.5 * diffusivity * runge_kutta_coefficients->implicit * dt,
        1.,
        buffers + i * NY,
        transposed + i * NY
    );
  }
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      scalar[j][i] += transposed[(i - 1) * NY + (j - 1)];
    }
  }
  if (0 != impose_boundary_condition_scalar_x(scalar_array)) {
    goto abort;
  }
  if (0 != impose_boundary_condition_scalar_y(scalar_array)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

