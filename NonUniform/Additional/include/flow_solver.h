#if !defined(FLOW_SOLVER_H)
#define FLOW_SOLVER_H

#include "array.h" // array_t
#include "domain.h" // domain_t
#include "tridiagonal_solver.h" // tridiagonal_t

typedef struct {
  // buffers to store intermediate data
  double * buf0;
  double * buf1;
  // x direction: gemm-related things
  double * forward_operator;
  double * backward_operator;
  double * eigenvalues;
  // y direction: tridiagonal_solver-related things
  tridiagonal_t * diagonals;
  double * tridiagonal_solver_buffers;
} poisson_solver_t;

typedef struct {
  // auxiliary buffers
  double * buffers;
  double * transposed;
  // matrices
  tridiagonal_t * x_diagonals;
  tridiagonal_t * y_diagonals;
} helmholtz_solver_t;

typedef struct {
  array_t * psi;
  array_t * dux;
  array_t * duy;
  array_t * dsalinity;
  array_t * dtemperature;
  array_t * dux_explicit_previous;
  array_t * duy_explicit_previous;
  array_t * dsalinity_explicit_previous;
  array_t * dtemperature_explicit_previous;
  poisson_solver_t poisson_solver;
  helmholtz_solver_t helmholtz_solver_ux;
  helmholtz_solver_t helmholtz_solver_uy;
  helmholtz_solver_t helmholtz_solver_scalar;
} flow_solver_t;

extern int flow_solver_init(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
);

extern int flow_solver_finalize(
    flow_solver_t * const flow_solver
);

#endif // FLOW_SOLVER_H
