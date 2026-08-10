#include "logger.h"
#include "./correct.h"

static int correct_ux(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    const flow_solver_t * const flow_solver,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  double * const * const psi = flow_solver->psi;
  double ** const ux = flow_field->ux;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      ux[j][i] -= dt / x_face_scale_factors[i] * (
          - psi[j    ][i - 1]
          + psi[j    ][i    ]
      );
    }
  }
  return 0;
}

static int correct_uy(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    const flow_solver_t * const flow_solver,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  double * const * const psi = flow_solver->psi;
  double ** const uy = flow_field->uy;
#pragma omp parallel for
  for (size_t j = 2; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      uy[j][i] -= dt / y_face_scale_factors[j] * (
          - psi[j - 1][i    ]
          + psi[j    ][i    ]
      );
    }
  }
  return 0;
}

int correct(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    const flow_solver_t * const flow_solver,
    const double dt
) {
  if (0 != correct_ux(domain, flow_field, flow_solver, dt)) {
    LOGGER_FAILURE("failed to correct ux");
    goto abort;
  }
  if (0 != correct_uy(domain, flow_field, flow_solver, dt)) {
    LOGGER_FAILURE("failed to correct uy");
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE("failed to correct velocity field");
  return 1;
}

