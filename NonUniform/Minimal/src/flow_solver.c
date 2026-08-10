#include <math.h>
#include "array.h"
#include "flow_solver.h"
#include "logger.h"
#include "memory.h"
#include "tridiagonal_solver.h"
#include "./eigendecomposition.h"

static int init_x_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const size_t nx = domain->x_coordinate.nitems;
  double ** const eigenvalues = &poisson_solver->eigenvalues;
  double ** const forward_operator = &poisson_solver->forward_operator;
  double ** const backward_operator = &poisson_solver->backward_operator;
  *eigenvalues = memory_alloc(nx, sizeof(double));
  *forward_operator = memory_alloc(nx * nx, sizeof(double));
  *backward_operator = memory_alloc(nx * nx, sizeof(double));
  // apply symmetrization and pass main/sub diagonal components to the eigendecomposition engine
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const double * const face_scale_factors = x_coordinate->face_scale_factors;
  const double * const center_scale_factors = x_coordinate->center_scale_factors;
  double * const main_diagonal = memory_alloc(nx, sizeof(double));
  double * const sub_diagonal = memory_alloc(nx - 1, sizeof(double));
  double * const eigenvectors_transposed = memory_alloc(nx * nx, sizeof(double));
  {
    const size_t n = 1;
    const double upper = 1. / center_scale_factors[n] / face_scale_factors[n + 1];
    main_diagonal[n - 1] = - upper;
  }
  for (size_t n = 2; n < nx; n++) {
    const double lower = 1. / center_scale_factors[n] / face_scale_factors[n    ];
    const double upper = 1. / center_scale_factors[n] / face_scale_factors[n + 1];
    main_diagonal[n - 1] = - lower - upper;
  }
  {
    const size_t n = nx;
    const double lower = 1. / center_scale_factors[n] / face_scale_factors[n    ];
    main_diagonal[n - 1] = - lower;
  }
  for (size_t n = 2; n < nx + 1; n++) {
    sub_diagonal[n - 2] = 1.
      / sqrt(center_scale_factors[n - 1])
      / sqrt(center_scale_factors[n])
      / face_scale_factors[n];
  }
  // find eigenvalues / eigenvectors
  if (0 != compute_eigendecomposition(
      nx * 8,
      nx,
      main_diagonal,
      sub_diagonal,
      *eigenvalues,
      eigenvectors_transposed
  )) {
    return 1;
  }
  // apply symmetrization operator to find forward / backward operators (between physical and spectral domains)
  // NOTE: notice the index difference (scale factors are defined on the domain whereas operators are pure matrices not recognizing the domain)
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < nx; j++) {
      (*forward_operator)[i * nx + j] = eigenvectors_transposed[i * nx + j] * sqrt(center_scale_factors[j + 1]);
    }
  }
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < nx; j++) {
      (*backward_operator)[i * nx + j] = 1. / sqrt(center_scale_factors[i + 1]) * eigenvectors_transposed[j * nx + i];
    }
  }
  // clean-up
  memory_free(main_diagonal);
  memory_free(sub_diagonal);
  memory_free(eigenvectors_transposed);
  return 0;
}

static int init_y_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
  double ** const tridiagonal_solver_lower_diagonals = &poisson_solver->tridiagonal_solver_lower_diagonals;
  double ** const tridiagonal_solver_main_diagonals = &poisson_solver->tridiagonal_solver_main_diagonals;
  double ** const tridiagonal_solver_upper_diagonals = &poisson_solver->tridiagonal_solver_upper_diagonals;
  double ** const tridiagonal_solver_buffers = &poisson_solver->tridiagonal_solver_buffers;
  *tridiagonal_solver_lower_diagonals = memory_alloc(ny, sizeof(double));
  *tridiagonal_solver_main_diagonals = memory_alloc(ny, sizeof(double));
  *tridiagonal_solver_upper_diagonals = memory_alloc(ny, sizeof(double));
  *tridiagonal_solver_buffers = memory_alloc(nx * ny, sizeof(double));
  for (size_t j = 1; j <= ny; j++) {
    double * const lower_diagonal = *tridiagonal_solver_lower_diagonals + j - 1;
    double * const main_diagonal = *tridiagonal_solver_main_diagonals + j - 1;
    double * const upper_diagonal = *tridiagonal_solver_upper_diagonals + j - 1;
    *lower_diagonal = 1. / y_center_scale_factors[j] / y_face_scale_factors[j    ];
    *upper_diagonal = 1. / y_center_scale_factors[j] / y_face_scale_factors[j + 1];
    *main_diagonal =
      1 == j ? - *upper_diagonal
      : ny == j ? - *lower_diagonal
      : - *lower_diagonal - *upper_diagonal;
  }
  return 0;
}

int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  // auxiliary buffers
  array_init(nx + 2, ny + 2, &flow_solver->psi);
  array_init(nx + 2, ny + 2, &flow_solver->dux);
  array_init(nx + 2, ny + 2, &flow_solver->duy);
  // poisson solver
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  double ** const buf0 = &poisson_solver->buf0;
  double ** const buf1 = &poisson_solver->buf1;
  *buf0 = memory_alloc(nx * ny, sizeof(double));
  *buf1 = memory_alloc(nx * ny, sizeof(double));
  // x direction: gemm-related things
  if (0 != init_x_solver(domain, poisson_solver)) {
    LOGGER_FAILURE("failed to initialise gemm part of poisson solver");
    goto abort;
  }
  // y direction: tridiagonal_solver-related things
  if (0 != init_y_solver(domain, poisson_solver)) {
    LOGGER_FAILURE("failed to initialise tridiagonal_solver part of poisson solver");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to initialise flow solver");
  return 1;
}

int flow_solver_finalize(
    flow_solver_t * const flow_solver
) {
  // auxiliary buffers
  array_finalize(&flow_solver->psi);
  array_finalize(&flow_solver->dux);
  array_finalize(&flow_solver->duy);
  // poisson solver
  poisson_solver_t * const poisson_solver = &flow_solver->poisson_solver;
  memory_free(poisson_solver->buf0);
  memory_free(poisson_solver->buf1);
  memory_free(poisson_solver->forward_operator);
  memory_free(poisson_solver->backward_operator);
  memory_free(poisson_solver->eigenvalues);
  memory_free(poisson_solver->tridiagonal_solver_lower_diagonals);
  memory_free(poisson_solver->tridiagonal_solver_main_diagonals);
  memory_free(poisson_solver->tridiagonal_solver_upper_diagonals);
  memory_free(poisson_solver->tridiagonal_solver_buffers);
  return 0;
}

