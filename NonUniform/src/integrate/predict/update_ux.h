#if !defined(UPDATE_UX_H)
#define UPDATE_UX_H

#include "array.h"
#include "flow_solver.h"
#include "param.h"
#include "runge_kutta.h"

extern int update_ux(
    const parameters_t * const parameters,
    const helmholtz_solver_t * const solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients,
    array_t * const dux_array,
    array_t * const ux_array
);

#endif // UPDATE_UX_H
