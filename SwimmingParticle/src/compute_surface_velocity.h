#if !defined(COMPUTE_SURFACE_VELOCITY_H)
#define COMPUTE_SURFACE_VELOCITY_H

#include "domain.h"

extern int compute_surface_velocity(
    const domain_t * const domain,
    const double * const surface_concentration,
    double * const surface_velocity
);

#endif // COMPUTE_SURFACE_VELOCITY_H
