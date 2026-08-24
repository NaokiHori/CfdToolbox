#include "boundary_condition.h"
#include "logger.h"
#include "tridiagonal_solver.h"
#include "./update_ux.h"

int update_ux(
    const parameters_t * const parameters,
    const helmholtz_solver_t * const solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients,
    array_t * const dux_array,
    array_t * const ux_array
) {
  const double diffusivity = parameters->fluid_diffusivity;
  double * const buffers = solver->buffers;
  double * const transposed = solver->transposed;
  const tridiagonal_t * const x_diagonals = solver->x_diagonals;
  const tridiagonal_t * const y_diagonals = solver->y_diagonals;
  double (* const dux)[NX + 2] = dux_array->buffer;
  double (* const ux)[NX + 2] = ux_array->buffer;
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    tridiagonal_solver_solve(
        NX - 1,
        x_diagonals,
        - 0.5 * diffusivity * runge_kutta_coefficients->implicit * dt,
        1.,
        buffers + (j - 1) * (NX - 1),
        &dux[j][2]
    );
  }
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 2; i <= NX; i++) {
      transposed[(i - 2) * NY + (j - 1)] = dux[j][i];
    }
  }
#pragma omp parallel for
  for (size_t i = 0; i < NX - 1; i++) {
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
    for (size_t i = 2; i <= NX; i++) {
      ux[j][i] += transposed[(i - 2) * NY + (j - 1)];
    }
  }
  if (0 != impose_boundary_condition_ux_x(ux_array)) {
    goto abort;
  }
  if (0 != impose_boundary_condition_ux_y(ux_array)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

