#include "boundary_condition.h"
#include "logger.h"
#include "tridiagonal_solver.h"
#include "./update_uy.h"

int update_uy(
    const parameters_t * const parameters,
    const helmholtz_solver_t * const solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients,
    array_t * const duy_array,
    array_t * const uy_array
) {
  const double diffusivity = parameters->fluid_diffusivity;
  double * const buffers = solver->buffers;
  double * const transposed = solver->transposed;
  const tridiagonal_t * const x_diagonals = solver->x_diagonals;
  const tridiagonal_t * const y_diagonals = solver->y_diagonals;
  double (* const duy)[NX + 2] = duy_array->buffer;
  double (* const uy)[NX + 2] = uy_array->buffer;
#pragma omp parallel for
  for (size_t j = 2; j <= NY; j++) {
    tridiagonal_solver_solve(
        NX,
        x_diagonals,
        - 0.5 * diffusivity * runge_kutta_coefficients->implicit * dt,
        1.,
        buffers + (j - 2) * NX,
        &duy[j][1]
    );
  }
#pragma omp parallel for
  for (size_t j = 2; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      transposed[(i - 1) * (NY - 1) + (j - 2)] = duy[j][i];
    }
  }
#pragma omp parallel for
  for (size_t i = 0; i < NX; i++) {
    tridiagonal_solver_solve(
        NY - 1,
        y_diagonals,
        - 0.5 * diffusivity * runge_kutta_coefficients->implicit * dt,
        1.,
        buffers + i * (NY - 1),
        transposed + i * (NY - 1)
    );
  }
#pragma omp parallel for
  for (size_t j = 2; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      uy[j][i] += transposed[(i - 1) * (NY - 1) + (j - 2)];
    }
  }
  if (0 != impose_boundary_condition_uy_x(uy_array)) {
    goto abort;
  }
  if (0 != impose_boundary_condition_uy_y(uy_array)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

