#if !defined(MONITOR_H)
#define MONITOR_H

#include "domain.h"

extern int monitor(
    const double l,
    const double kappa,
    const vector_t * const velocity,
    const scalar_t * const temperature,
    double * const max_velocity,
    double * const nusselt
);

#endif // MONITOR_H
