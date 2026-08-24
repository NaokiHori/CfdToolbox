#include <stddef.h>
#include "logger.h"
#include "memory.h"
#include "./helmholtz_solver_ux.h"

static int init_x_diagonals(
    const coordinate_t * const coordinate,
    tridiagonal_t ** const diagonals
) {
  const tridiagonal_t * const face_laplace_operators = coordinate->face_laplace_operators;
  *diagonals = memory_alloc(NX - 1, sizeof(tridiagonal_t));
  for (size_t i = 2; i <= NX; i++) {
    tridiagonal_t * const diagonal = *diagonals + i - 2;
    *diagonal = face_laplace_operators[i];
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

int init_helmholtz_solver_ux(
    const domain_t * const domain,
    helmholtz_solver_t * const solver
) {
  solver->buffers = memory_alloc((NX - 1) * NY, sizeof(double));
  solver->transposed = memory_alloc((NX - 1) * NY, sizeof(double));
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

