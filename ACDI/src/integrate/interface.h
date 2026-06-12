#if !defined(INTERFACE_H)
#define INTERFACE_H

#include "domain.h"
#include "flow_field.h"
#include "interface_field.h"
#include "interface_solver.h"

extern int update_interface_field(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    interface_field_t * const interface_field,
    interface_solver_t * const interface_solver,
    const double dt
);

#endif // INTERFACE_H
