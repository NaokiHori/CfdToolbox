#if !defined(COMPUTE_DUY_H)
#define COMPUTE_DUY_H

#include "domain.h"
#include "flow_field.h"
#include "interface_field.h"
#include "interface_solver.h"

int compute_duy(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const interface_field_t * const interface_field,
    const interface_solver_t * const interface_solver,
    const double dt,
    double ** const duy
);

#endif // COMPUTE_DUY_H
