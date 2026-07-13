#if !defined(IMMERSED_BOUNDARY_H)
#define IMMERSED_BOUNDARY_H

#include <stdbool.h> // bool
#include "domain.h" // domain_t
#include "flow_field.h" // flow_field_t
#include "particle_stuffs.h" // coordinate_t, particle_stuffs_t

extern int immersed_boundary_reset_particle_increments(
    particle_stuffs_t * const particle_stuffs
);

extern int immersed_boundary_compute_internal_momentum(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    particle_stuffs_t * const particle_stuffs,
    const bool update_new
);

extern int immersed_boundary_exchange_momentum(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    particle_stuffs_t * const particle_stuffs,
    const double dt
);

extern int immersed_boundary_correct_velocity(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    particle_stuffs_t * const particle_stuffs,
    const double dt
);

extern int immersed_boundary_update_particle_increments(
    particle_stuffs_t * const particle_stuffs,
    const double dt,
    double * const tolerance
);

extern int immersed_boundary_update_particle(
    const domain_t * const domain,
    particle_stuffs_t * const particle_stuffs
);

#endif // IMMERSED_BOUNDARY_H
