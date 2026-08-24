#include <math.h>
#include "array.h"
#include "flow_solver.h"
#include "logger.h"
#include "memory.h"
#include "param.h"
#include "tridiagonal_solver.h"
#include "./flow_solver/helmholtz_solver_scalar.h"
#include "./flow_solver/helmholtz_solver_ux.h"
#include "./flow_solver/helmholtz_solver_uy.h"
#include "./flow_solver/poisson_solver.h"

int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
) {
  // auxiliary buffers
  flow_solver->psi = memory_alloc(1, sizeof(array_t));
  flow_solver->dux = memory_alloc(1, sizeof(array_t));
  flow_solver->duy = memory_alloc(1, sizeof(array_t));
  flow_solver->dtemperature = memory_alloc(1, sizeof(array_t));
  flow_solver->dsalinity = memory_alloc(1, sizeof(array_t));
  flow_solver->dux_explicit_previous = memory_alloc(1, sizeof(array_t));
  flow_solver->duy_explicit_previous = memory_alloc(1, sizeof(array_t));
  flow_solver->dtemperature_explicit_previous = memory_alloc(1, sizeof(array_t));
  flow_solver->dsalinity_explicit_previous = memory_alloc(1, sizeof(array_t));
  // poisson solver
  if (0 != init_poisson_solver(domain, &flow_solver->poisson_solver)) {
    goto abort;
  }
  // helmholtz solvers
  if (0 != init_helmholtz_solver_ux(domain, &flow_solver->helmholtz_solver_ux)) {
    goto abort;
  }
  if (0 != init_helmholtz_solver_uy(domain, &flow_solver->helmholtz_solver_uy)) {
    goto abort;
  }
  if (0 != init_helmholtz_solver_scalar(domain, &flow_solver->helmholtz_solver_scalar)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

int flow_solver_finalize(
    flow_solver_t * const flow_solver
) {
  // auxiliary buffers
  memory_free(flow_solver->psi);
  memory_free(flow_solver->dux);
  memory_free(flow_solver->duy);
  memory_free(flow_solver->dtemperature);
  memory_free(flow_solver->dsalinity);
  memory_free(flow_solver->dux_explicit_previous);
  memory_free(flow_solver->duy_explicit_previous);
  memory_free(flow_solver->dtemperature_explicit_previous);
  memory_free(flow_solver->dsalinity_explicit_previous);
  // poisson solver
  {
    poisson_solver_t * const solver = &flow_solver->poisson_solver;
    memory_free(solver->buf0);
    memory_free(solver->buf1);
    memory_free(solver->forward_operator);
    memory_free(solver->backward_operator);
    memory_free(solver->eigenvalues);
    memory_free(solver->tridiagonal_solver_buffers);
    memory_free(solver->diagonals);
  }
  // helmholtz solvers
  {
    helmholtz_solver_t * const solver = &flow_solver->helmholtz_solver_ux;
    memory_free(solver->x_diagonals);
    memory_free(solver->y_diagonals);
    memory_free(solver->buffers);
    memory_free(solver->transposed);
  }
  {
    helmholtz_solver_t * const solver = &flow_solver->helmholtz_solver_uy;
    memory_free(solver->x_diagonals);
    memory_free(solver->y_diagonals);
    memory_free(solver->buffers);
    memory_free(solver->transposed);
  }
  {
    helmholtz_solver_t * const solver = &flow_solver->helmholtz_solver_scalar;
    memory_free(solver->x_diagonals);
    memory_free(solver->y_diagonals);
    memory_free(solver->buffers);
    memory_free(solver->transposed);
  }
  return 0;
}

