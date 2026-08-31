#if !defined(DECIDE_DT_H)
#define DECIDE_DT_H

#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "param.h"

extern int decide_dt(
    const parameters_t * const parameters,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
);

#endif // DECIDE_DT_H
