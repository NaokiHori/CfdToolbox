#if !defined(UPDATE_UY_H)
#define UPDATE_UY_H

#include "array.h"
#include "flow_solver.h"
#include "param.h"
#include "runge_kutta.h"

extern int update_uy(
    const parameters_t * const parameters,
    const helmholtz_solver_t * const solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients,
    array_t * const duy_array,
    array_t * const uy_array
);

#endif // UPDATE_UY_H
