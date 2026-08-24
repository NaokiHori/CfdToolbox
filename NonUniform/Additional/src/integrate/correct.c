#include "logger.h"
#include "./correct.h"

static int correct_ux(
    const domain_t * const domain,
    const double dt,
    const array_t * const psi_array,
    array_t * const ux_array
) {
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  const double (* const psi)[NX + 2] = psi_array->buffer;
  double (* const ux)[NX + 2] = ux_array->buffer;
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 2; i <= NX; i++) {
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
    const double dt,
    const array_t * const psi_array,
    array_t * const uy_array
) {
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double (* const psi)[NX + 2] = psi_array->buffer;
  double (* const uy)[NX + 2] = uy_array->buffer;
#pragma omp parallel for
  for (size_t j = 2; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
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
    const double dt,
    const flow_solver_t * const flow_solver,
    flow_field_t * const flow_field
) {
  const array_t * const psi_array = flow_solver->psi;
  array_t * const ux_array = flow_field->ux;
  array_t * const uy_array = flow_field->uy;
  if (0 != correct_ux(domain, dt, psi_array, ux_array)) {
    goto abort;
  }
  if (0 != correct_uy(domain, dt, psi_array, uy_array)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

