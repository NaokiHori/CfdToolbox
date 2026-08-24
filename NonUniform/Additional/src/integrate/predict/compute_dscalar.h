#if !defined(COMPUTE_DSCALAR_H)
#define COMPUTE_DSCALAR_H

#include "array.h"
#include "domain.h"
#include "flow_field.h"
#include "runge_kutta.h"

int compute_dscalar(
    const double diffusivity,
    const domain_t * const domain,
    const array_t * const ux_array,
    const array_t * const uy_array,
    const array_t * const scalar_array,
    array_t * const dscalar,
    array_t * const dscalar_explicit_previous,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients
);

#endif // COMPUTE_DSCALAR_H
