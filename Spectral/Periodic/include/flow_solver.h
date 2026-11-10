#if !defined(FLOW_SOLVER_H)
#define FLOW_SOLVER_H

#include <complex.h>
#include "domain.h"
#include "flow_field.h"
#include "runge_kutta.h"
#include "transform.h"

#define N_BUFFERS 4

typedef struct {
  double * integrating_factor_coefficients;
  // runge-kutta integrators
  runge_kutta_t * runge_kutta_ux;
  runge_kutta_t * runge_kutta_uy;
  // auxiliary buffers
  void * buffers[N_BUFFERS];
} flow_solver_t;

extern int initialize_flow_solver(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
);

extern int finalize_flow_solver(
    flow_solver_t * const flow_solver
);

extern int integrate(
    const domain_t * const domain,
    transformer_t * const transformer,
    flow_solver_t * const flow_solver,
    flow_field_t * const flow_field,
    double * const dt
);

#endif // FLOW_SOLVER_H
