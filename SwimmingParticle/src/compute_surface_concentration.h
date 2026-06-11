#if !defined(COMPUTE_SURFACE_CONCENTRATION_H)
#define COMPUTE_SURFACE_CONCENTRATION_H

#include "domain.h"
#include "rdft.h"

extern int compute_surface_concentration(
    const domain_t * const domain,
    rdft_plan_t * const rdft_plan,
    double ** const concentration,
    double * const surface_concentration
);

#endif // COMPUTE_SURFACE_CONCENTRATION_H
