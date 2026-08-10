#include "./update_pressure.h"

int update_pressure(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  double ** const psi = flow_solver->psi;
  double ** const p = flow_field->p;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      p[j][i] += psi[j][i];
    }
  }
  return 0;
}

