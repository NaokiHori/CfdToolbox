#include "array.h"
#include "flow_solver.h"
#include "load_grid.h"
#include "logger.h"
#include "memory.h"
#include "tridiagonal_solver.h"

static int init_x_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const size_t nx = domain->nx;
  double ** const eigenvalues = &poisson_solver->eigenvalues;
  double ** const forward_operator = &poisson_solver->forward_operator;
  double ** const backward_operator = &poisson_solver->backward_operator;
  *eigenvalues = memory_alloc(nx, sizeof(double));
  *forward_operator = memory_alloc(nx * nx, sizeof(double));
  *backward_operator = memory_alloc(nx * nx, sizeof(double));
  if (0 != try_load_npy(
      GRID_DIRECTORY "x_eigenvalues.npy",
      1,
      (size_t [1]){nx},
      "'<f8'",
      false,
      sizeof(double),
      *eigenvalues
  )) {
    return 1;
  }
  if (0 != try_load_npy(
      GRID_DIRECTORY "x_forward_operator.npy",
      2,
      (size_t [2]){nx, nx},
      "'<f8'",
      false,
      sizeof(double),
      *forward_operator
  )) {
    return 1;
  }
  if (0 != try_load_npy(
      GRID_DIRECTORY "x_backward_operator.npy",
      2,
      (size_t [2]){nx, nx},
      "'<f8'",
      false,
      sizeof(double),
      *backward_operator
  )) {
    return 1;
  }
  return 0;
}

static int init_y_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double * const y_center_scale_factors = domain->y_coordinate.center_scale_factors;
  tridiagonal_solver_plan_t ** const tridiagonal_solver_plan = &poisson_solver->tridiagonal_solver_plan;
  if (0 != tridiagonal_solver_init_plan(ny, nx, tridiagonal_solver_plan)) {
    LOGGER_FAILURE("failed to initialise tridiagonal_solver solver");
    goto abort;
  }
  double ** const tridiagonal_solver_l = &poisson_solver->tridiagonal_solver_l;
  double ** const tridiagonal_solver_m = &poisson_solver->tridiagonal_solver_m;
  double ** const tridiagonal_solver_u = &poisson_solver->tridiagonal_solver_u;
  *tridiagonal_solver_l = memory_alloc(ny, sizeof(double));
  *tridiagonal_solver_m = memory_alloc(ny, sizeof(double));
  *tridiagonal_solver_u = memory_alloc(ny, sizeof(double));
  for (size_t j = 1; j <= ny; j++) {
    const double l = 1. / y_center_scale_factors[j] / y_face_scale_factors[j    ];
    const double u = 1. / y_center_scale_factors[j] / y_face_scale_factors[j + 1];
    (*tridiagonal_solver_l)[j - 1] = + 1. * l;
    (*tridiagonal_solver_u)[j - 1] = + 1. * u;
    (*tridiagonal_solver_m)[j - 1] =
      1 == j ? - 1. * u
      : ny == j ? - 1. * l
      : - 1. * l - 1. * u;
  }
  return 0;
abort:
  return 1;
}

int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
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
  // x direction: dft-related things
  if (0 != init_x_solver(domain, poisson_solver)) {
    LOGGER_FAILURE("failed to initialise dft part of poisson solver");
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
  tridiagonal_solver_destroy_plan(&poisson_solver->tridiagonal_solver_plan);
  memory_free(poisson_solver->buf0);
  memory_free(poisson_solver->buf1);
  memory_free(poisson_solver->forward_operator);
  memory_free(poisson_solver->backward_operator);
  memory_free(poisson_solver->eigenvalues);
  memory_free(poisson_solver->tridiagonal_solver_l);
  memory_free(poisson_solver->tridiagonal_solver_m);
  memory_free(poisson_solver->tridiagonal_solver_u);
  return 0;
}

