#include "boundary_condition.h"
#include "logger.h"
#include "memory.h"
#include "tridiagonal_solver.h"
#include "./predict.h"
#include "./predict/compute_dscalar.h"
#include "./predict/compute_dux.h"
#include "./predict/compute_duy.h"
#include "./predict/update_scalar.h"
#include "./predict/update_ux.h"
#include "./predict/update_uy.h"

int predict(
    const parameters_t * const parameters,
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients
) {
  // compute increments
  if (0 != compute_dux(
      parameters,
      domain,
      flow_field,
      flow_solver,
      dt,
      runge_kutta_coefficients
  )) {
    goto abort;
  }
  if (0 != compute_duy(
      parameters,
      domain,
      flow_field,
      flow_solver,
      dt,
      runge_kutta_coefficients
  )) {
    goto abort;
  }
  if (0 != compute_dscalar(
      parameters->temperature_diffusivity,
      domain,
      flow_field->ux,
      flow_field->uy,
      flow_field->temperature,
      flow_solver->dtemperature,
      flow_solver->dtemperature_explicit_previous,
      dt,
      runge_kutta_coefficients
  )) {
    goto abort;
  }
  if (0 != compute_dscalar(
      parameters->salinity_diffusivity,
      domain,
      flow_field->ux,
      flow_field->uy,
      flow_field->salinity,
      flow_solver->dsalinity,
      flow_solver->dsalinity_explicit_previous,
      dt,
      runge_kutta_coefficients
  )) {
    goto abort;
  }
  // update field
  if (0 != update_ux(
      parameters,
      &flow_solver->helmholtz_solver_ux,
      dt,
      runge_kutta_coefficients,
      flow_solver->dux,
      flow_field->ux
  )) {
    goto abort;
  }
  if (0 != update_uy(
      parameters,
      &flow_solver->helmholtz_solver_uy,
      dt,
      runge_kutta_coefficients,
      flow_solver->duy,
      flow_field->uy
  )) {
    goto abort;
  }
  if (0 != update_scalar(
      parameters->temperature_diffusivity,
      &flow_solver->helmholtz_solver_scalar,
      dt,
      runge_kutta_coefficients,
      flow_solver->dtemperature,
      flow_field->temperature
  )) {
    goto abort;
  }
  if (0 != update_scalar(
      parameters->salinity_diffusivity,
      &flow_solver->helmholtz_solver_scalar,
      dt,
      runge_kutta_coefficients,
      flow_solver->dsalinity,
      flow_field->salinity
  )) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

