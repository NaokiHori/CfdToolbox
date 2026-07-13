#include <math.h>
#include "memory.h"
#include "array.h"
#include "logger.h"
#include "interface_solver.h"

int interface_solver_init(
    const domain_t * const domain,
    interface_solver_t * const interface_solver
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  array_init(nx + 2, ny + 2, &interface_solver->sdf);
  array_init(nx + 2, ny + 2, &interface_solver->flux_x);
  array_init(nx + 2, ny + 2, &interface_solver->flux_y);
  array_init(nx + 2, ny + 2, &interface_solver->normal_x);
  array_init(nx + 2, ny + 2, &interface_solver->normal_y);
  array_init(nx + 2, ny + 2, &interface_solver->curvature);
  return 0;
}

int interface_solver_finalize(
    interface_solver_t * const interface_solver
) {
  array_finalize(&interface_solver->sdf);
  array_finalize(&interface_solver->flux_x);
  array_finalize(&interface_solver->flux_y);
  array_finalize(&interface_solver->normal_x);
  array_finalize(&interface_solver->normal_y);
  array_finalize(&interface_solver->curvature);
  return 0;
}

