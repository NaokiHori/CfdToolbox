#if !defined(PREDICT_H)
#define PREDICT_H

#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "interface_field.h"
#include "interface_solver.h"

extern int predict(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    flow_solver_t * const flow_solver,
    interface_field_t * const interface_field,
    interface_solver_t * const interface_solver,
    const double dt
);

#endif // PREDICT_H
