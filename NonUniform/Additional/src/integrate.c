#include "logger.h"
#include "runge_kutta.h"
#include "./integrate.h"
#include "./integrate/correct.h"
#include "./integrate/decide_dt.h"
#include "./integrate/predict.h"
#include "./integrate/solve_poisson.h"
#include "./integrate/update_pressure.h"

int integrate(
    const parameters_t * const parameters,
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    double * const dt
) {
  if (0 != decide_dt(parameters, domain, flow_field, dt)) {
    return 1;
  }
  for (size_t runge_kutta_step = 0; runge_kutta_step < N_RUNGE_KUTTA_STEP; runge_kutta_step++) {
    runge_kutta_coefficients_t runge_kutta_coefficients = {0};
    if (0 != runge_kutta_coefficients_setup(runge_kutta_step, &runge_kutta_coefficients)) {
      goto abort;
    }
    if (0 != predict(
        parameters,
        domain,
        flow_field,
        flow_solver,
        *dt,
        &runge_kutta_coefficients
    )) {
      goto abort;
    }
    if (0 != solve_poisson(domain, flow_field, flow_solver, runge_kutta_coefficients.implicit * *dt)) {
      goto abort;
    }
    if (0 != correct(
        domain,
        runge_kutta_coefficients.implicit * *dt,
        flow_solver,
        flow_field
    )) {
      goto abort;
    }
    if (0 != update_pressure(
        parameters,
        domain,
        runge_kutta_coefficients.implicit * *dt,
        flow_solver,
        flow_field
    )) {
      goto abort;
    }
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

