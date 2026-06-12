#if !defined(INTERFACE_SOLVER_H)
#define INTERFACE_SOLVER_H

#include "domain.h" // domain_t

typedef struct {
  // signed distance function
  double ** sdf;
  // fluxes to update vof
  double ** flux_x;
  double ** flux_y;
  // normal vectors
  double ** normal_x;
  double ** normal_y;
  //
  double ** curvature;
} interface_solver_t;

extern int interface_solver_init(
    const domain_t * const domain,
    interface_solver_t * const interface_solver
);

extern int interface_solver_finalize(
    interface_solver_t * const interface_solver
);

#endif // INTERFACE_SOLVER_H
