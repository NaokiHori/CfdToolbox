#include <stddef.h>
#include "logger.h"
#include "memory.h"
#include "./helmholtz_solver_scalar.h"

static int init_x_diagonals(
    const coordinate_t * const coordinate,
    tridiagonal_t ** const diagonals
) {
  const tridiagonal_t * const center_laplace_operators = coordinate->center_laplace_operators;
  *diagonals = memory_alloc(NX, sizeof(tridiagonal_t));
  for (size_t i = 1; i <= NX; i++) {
    tridiagonal_t * const diagonal = *diagonals + i - 1;
    *diagonal = center_laplace_operators[i];
    double * const main = &diagonal->main;
    const double lower = diagonal->lower;
    const double upper = diagonal->upper;
    // impose neumann boundary conditions
    if (1 == i) {
      *main += lower;
    }
    if (NX == i) {
      *main += upper;
    }
  }
  return 0;
}

static int init_y_diagonals(
    const coordinate_t * const coordinate,
    tridiagonal_t ** const diagonals
) {
  const tridiagonal_t * const center_laplace_operators = coordinate->center_laplace_operators;
  *diagonals = memory_alloc(NY, sizeof(tridiagonal_t));
  for (size_t j = 1; j <= NY; j++) {
    tridiagonal_t * const diagonal = *diagonals + j - 1;
    *diagonal = center_laplace_operators[j];
  }
  return 0;
}

int init_helmholtz_solver_scalar(
    const domain_t * const domain,
    helmholtz_solver_t * const solver
) {
  solver->buffers = memory_alloc(NX * NY, sizeof(double));
  solver->transposed = memory_alloc(NX * NY, sizeof(double));
  if (0 != init_x_diagonals(&domain->x_coordinate, &solver->x_diagonals)) {
    goto abort;
  }
  if (0 != init_y_diagonals(&domain->y_coordinate, &solver->y_diagonals)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

