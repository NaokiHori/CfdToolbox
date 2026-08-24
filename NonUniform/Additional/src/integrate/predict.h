#if !defined(PREDICT_H)
#define PREDICT_H

#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "param.h"
#include "runge_kutta.h"

extern int predict(
    const parameters_t * const parameters,
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    const double dt,
    const runge_kutta_coefficients_t * const runge_kutta_coefficients
);

#endif // PREDICT_H
