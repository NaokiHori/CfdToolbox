#if !defined(PARTICLE_H)
#define PARTICLE_H

#include <stddef.h>
#include "domain.h"

typedef struct {
  double current_velocity[NDIMS];
  double updated_velocity[NDIMS];
  double position[NDIMS];
} particle_t;

extern int compute_particle_velocity(
    const domain_t * const domain,
    double * const surface_concentration,
    double * velocity
);

extern int update_particle_position(
    const double dt,
    const size_t runge_kutta_step,
    particle_t * const particle
);

#endif // PARTICLE_H
