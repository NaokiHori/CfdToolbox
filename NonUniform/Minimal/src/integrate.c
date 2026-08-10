#include "./integrate.h"
#include "./integrate/correct.h"
#include "./integrate/decide_dt.h"
#include "./integrate/predict.h"
#include "./integrate/solve_poisson.h"
#include "./integrate/update_pressure.h"

int integrate(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    double * const dt
) {
  if (0 != decide_dt(domain, flow_field, dt)) {
    return 1;
  }
  if (0 != predict(domain, flow_field, flow_solver, *dt)) {
    return 1;
  }
  if (0 != solve_poisson(domain, flow_field, flow_solver, *dt)) {
    return 1;
  }
  if (0 != correct(domain, flow_field, flow_solver, *dt)) {
    return 1;
  }
  if (0 != update_pressure(domain, flow_field, flow_solver)) {
    return 1;
  }
  return 0;
}

