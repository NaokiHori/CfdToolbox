#include <math.h>
#include "logger.h"
#include "./immersed_boundary.h"
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
    particle_stuffs_t * const particle_stuffs,
    double * const dt
) {
  if (0 != decide_dt(domain, flow_field, dt)) {
    return 1;
  }
  if (0 != immersed_boundary_reset_particle_increments(particle_stuffs)) {
    return 1;
  }
  if (0 != immersed_boundary_compute_internal_momentum(
      domain,
      flow_field,
      particle_stuffs,
      false
  )) {
    return 1;
  }
  if (0 != predict(domain, flow_field, flow_solver, *dt)) {
    return 1;
  }
  if (0 != immersed_boundary_exchange_momentum(
      domain,
      flow_field,
      particle_stuffs,
      *dt
  )) {
    return 1;
  }
  if (0 != immersed_boundary_correct_velocity(
      domain,
      flow_field,
      particle_stuffs,
      *dt
  )) {
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
  double tolerance = 0.;
  for (;;) {
    if (0 != immersed_boundary_compute_internal_momentum(
        domain,
        flow_field,
        particle_stuffs,
        true
    )) {
      return 1;
    }
    if (0 != immersed_boundary_update_particle_increments(
        particle_stuffs,
        *dt,
        &tolerance
    )) {
      return 1;
    }
    if (tolerance < particle_stuffs->tolerance_threshold) {
      break;
    }
  }
  if (0 != immersed_boundary_update_particle(domain, particle_stuffs)) {
    return 1;
  }
  return 0;
}

