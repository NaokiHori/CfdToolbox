#if !defined(UPDATE_SCALAR_H)
#define UPDATE_SCALAR_H

#include "array.h"
#include "flow_solver.h"
#include "runge_kutta.h"

extern int update_scalar(
    const double diffusivity,
    const helmholtz_solver_t * const solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients,
    array_t * const dscalar_array,
    array_t * const scalar_array
);

#endif // UPDATE_SCALAR_H
