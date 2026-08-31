#include <stddef.h> // size_t
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "param.h"
#include "./integrate.h"
#include "./monitor.h"
#include "./save.h"

typedef struct {
  double monitor;
  double save;
} schedule_t;

int main(
    void
) {
  parameters_t parameters = {0};
  domain_t domain = {0};
  flow_field_t flow_field = {0};
  flow_solver_t flow_solver = {0};
  if (0 != setup_parameters(
      &(non_dimensional_numbers_t) {
        .density_ratio = 2e+0,
        .lewis = 5e+0,
        .prandtl = 4e+0,
        .rayleigh = 1e+7,
      },
      &parameters
  )) {
    return 1;
  }
  if (0 != domain_init(&domain)) {
    return 1;
  }
  if (0 != flow_field_init(&domain, &flow_field)) {
    return 1;
  }
  if (0 != flow_solver_init(&domain, &flow_solver)) {
    return 1;
  }
  const double time_max = 1.e+2;
  const schedule_t rate = {
    .monitor = 5.e-1,
    .save = 5.e+0,
  };
  schedule_t next = rate;
  for (double time = 0.; time < time_max; ) {
    static size_t step = 0;
    double dt = 0.;
    if (0 != integrate(&parameters, &domain, &flow_field, &flow_solver, &dt)) {
      break;
    }
    step += 1;
    time += dt;
    if (next.monitor < time) {
      monitor(step, time, dt, &domain, &flow_field);
      next.monitor += rate.monitor;
    }
    if (next.save < time) {
      static size_t id = 0;
      save(id, step, time, &domain, &flow_field);
      id += 1;
      next.save += rate.save;
    }
  }
  if (0 != domain_finalize(&domain)) {
    return 1;
  }
  if (0 != flow_field_finalize(&flow_field)) {
    return 1;
  }
  if (0 != flow_solver_finalize(&flow_solver)) {
    return 1;
  }
  return 0;
}

