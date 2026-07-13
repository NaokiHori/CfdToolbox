#if !defined(INTEGRATE_H)
#define INTEGRATE_H

#include "domain.h" // domain_t
#include "flow_field.h" // flow_field_t
#include "flow_solver.h" // flow_solver_t
#include "particle_stuffs.h" // particle_stuffs_t

extern int integrate(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    particle_stuffs_t * const particle_stuffs,
    double * const dt
);

#endif // INTEGRATE_H
