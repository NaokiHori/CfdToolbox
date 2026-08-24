#if !defined(HELMHOLTZ_SOLVER_SCALAR_H)
#define HELMHOLTZ_SOLVER_SCALAR_H

#include "domain.h"
#include "flow_solver.h"

extern int init_helmholtz_solver_scalar(
    const domain_t * const domain,
    helmholtz_solver_t * const solver
);

#endif // HELMHOLTZ_SOLVER_SCALAR_H
