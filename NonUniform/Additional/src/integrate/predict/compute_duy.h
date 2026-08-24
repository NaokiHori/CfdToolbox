#if !defined(COMPUTE_DUY_H)
#define COMPUTE_DUY_H

#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "param.h"
#include "runge_kutta.h"

int compute_duy(
    const parameters_t * const parameters,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const flow_solver_t * const flow_solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients
);

#endif // COMPUTE_DUY_H
