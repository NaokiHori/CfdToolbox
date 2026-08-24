#include <math.h>
#include "logger.h"
#include "memory.h"
#include "./eigendecomposition.h"
#include "./poisson_solver.h"

static int init_poisson_solver_x(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  double ** const eigenvalues = &poisson_solver->eigenvalues;
  double ** const forward_operator = &poisson_solver->forward_operator;
  double ** const backward_operator = &poisson_solver->backward_operator;
  *eigenvalues = memory_alloc(NX, sizeof(double));
  *forward_operator = memory_alloc(NX * NX, sizeof(double));
  *backward_operator = memory_alloc(NX * NX, sizeof(double));
  // apply symmetrization and pass main/sub diagonal components to the eigendecomposition engine
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const tridiagonal_t * const center_laplace_operators = x_coordinate->center_laplace_operators;
  const double * const face_scale_factors = x_coordinate->face_scale_factors;
  const double * const center_scale_factors = x_coordinate->center_scale_factors;
  double * const main_diagonals = memory_alloc(NX, sizeof(double));
  double * const sub_diagonals = memory_alloc(NX - 1, sizeof(double));
  double * const eigenvectors_transposed = memory_alloc(NX * NX, sizeof(double));
  for (size_t n = 1; n <= NX; n++) {
    double * const main_diagonal = main_diagonals + n - 1;
    const tridiagonal_t * const center_laplace_operator = center_laplace_operators + n;
    *main_diagonal = center_laplace_operator->main;
    // impose neumann boundary conditions
    if (1 == n) {
      *main_diagonal += center_laplace_operator->lower;
    }
    if (NX == n) {
      *main_diagonal += center_laplace_operator->upper;
    }
  }
  for (size_t n = 2; n < NX + 1; n++) {
    sub_diagonals[n - 2] = 1.
      / sqrt(center_scale_factors[n - 1])
      / sqrt(center_scale_factors[n])
      / face_scale_factors[n];
  }
  // find eigenvalues / eigenvectors
  if (0 != compute_eigendecomposition(
      8,
      NX,
      main_diagonals,
      sub_diagonals,
      *eigenvalues,
      eigenvectors_transposed
  )) {
    goto abort;
  }
  // apply symmetrization operator to find forward / backward operators (between physical and spectral domains)
  // NOTE: notice the index difference (scale factors are defined on the domain whereas operators are pure matrices not recognizing the domain)
  for (size_t i = 0; i < NX; i++) {
    for (size_t j = 0; j < NX; j++) {
      (*forward_operator)[i * NX + j] = eigenvectors_transposed[i * NX + j] * sqrt(center_scale_factors[j + 1]);
    }
  }
  for (size_t i = 0; i < NX; i++) {
    for (size_t j = 0; j < NX; j++) {
      (*backward_operator)[i * NX + j] = 1. / sqrt(center_scale_factors[i + 1]) * eigenvectors_transposed[j * NX + i];
    }
  }
  // clean-up
  memory_free(main_diagonals);
  memory_free(sub_diagonals);
  memory_free(eigenvectors_transposed);
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

static int init_poisson_solver_y(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const coordinate_t * const y_coordinate = &domain->y_coordinate;
  const tridiagonal_t * const center_laplace_operators = y_coordinate->center_laplace_operators;
  tridiagonal_t ** const diagonals = &poisson_solver->diagonals;
  double ** const tridiagonal_solver_buffers = &poisson_solver->tridiagonal_solver_buffers;
  *diagonals = memory_alloc(NY, sizeof(tridiagonal_t));
  *tridiagonal_solver_buffers = memory_alloc(NX * NY, sizeof(double));
  for (size_t j = 1; j <= NY; j++) {
    tridiagonal_t * const diagonal = *diagonals + j - 1;
    *diagonal = center_laplace_operators[j];
    double * const main = &diagonal->main;
    // impose neumann boundary conditions
    if (1 == j) {
      *main += diagonal->lower;
    }
    if (NY == j) {
      *main += diagonal->upper;
    }
  }
  return 0;
}

int init_poisson_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  double ** const buf0 = &poisson_solver->buf0;
  double ** const buf1 = &poisson_solver->buf1;
  *buf0 = memory_alloc(NX * NY, sizeof(double));
  *buf1 = memory_alloc(NX * NY, sizeof(double));
  // x direction: gemm-related things
  if (0 != init_poisson_solver_x(domain, poisson_solver)) {
    goto abort;
  }
  // y direction: tridiagonal-solver-related things
  if (0 != init_poisson_solver_y(domain, poisson_solver)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

