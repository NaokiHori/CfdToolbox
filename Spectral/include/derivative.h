#if !defined(DERIVATIVE_H)
#define DERIVATIVE_H

#include <complex.h>
#include "domain.h"

extern int compute_x_derivative(
    const domain_t * const domain,
    const double sign,
    const double complex * const src,
    double complex * const dst
);

extern int compute_y_derivative(
    const domain_t * const domain,
    const double sign,
    const double complex * const src,
    double complex * const dst
);

#endif // DERIVATIVE_H
