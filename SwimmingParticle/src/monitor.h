#if !defined(MONITOR_H)
#define MONITOR_H

#include "domain.h"
#include "particle.h"

extern int monitor(
    const domain_t * const domain,
    const double time,
    double ** const concentration,
    const particle_t * const particle
);

#endif // MONITOR_H
