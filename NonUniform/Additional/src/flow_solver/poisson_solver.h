#if !defined(POISSON_SOLVER_H)
#define POISSON_SOLVER_H

#include "domain.h"
#include "flow_solver.h"

extern int init_poisson_solver(
    const domain_t * const domain,
    poisson_solver_t * const poisson_solver
);

#endif // POISSON_SOLVER_H
