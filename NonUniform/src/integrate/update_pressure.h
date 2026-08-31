#if !defined(UPDATE_PRESSURE_H)
#define UPDATE_PRESSURE_H

#include "array.h"
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "param.h"

extern int update_pressure(
    const parameters_t * const parameters,
    const domain_t * const domain,
    const double dt,
    const flow_solver_t * const flow_solver,
    flow_field_t * const flow_field
);

#endif // UPDATE_PRESSURE_H
