#if !defined(FLOW_SOLVER_H)
#define FLOW_SOLVER_H

#include "domain.h" // domain_t
#include "tridiagonal_solver.h" // tridiagonal_solver_plan_t

// variables used to solve Poisson equations
typedef struct {
  // buffers to store intermediate data
  double * buf0;
  double * buf1;
  // x direction: gemm-related things
  double * forward_operator;
  double * backward_operator;
  double * eigenvalues;
  // y direction: tridiagonal_solver-related things
  tridiagonal_solver_plan_t * tridiagonal_solver_plan;
  double * tridiagonal_solver_l;
  double * tridiagonal_solver_m;
  double * tridiagonal_solver_u;
} poisson_solver_t;

typedef struct {
  double ** psi;
  double ** dux;
  double ** duy;
  poisson_solver_t poisson_solver;
} flow_solver_t;

extern int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
);

extern int flow_solver_finalize(
    flow_solver_t * const flow_solver
);

#endif // FLOW_SOLVER_H
