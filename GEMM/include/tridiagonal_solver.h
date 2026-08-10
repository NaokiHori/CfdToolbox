#if !defined(TRIDIAGONAL_SOLVER_H)
#define TRIDIAGONAL_SOLVER_H

#include <stddef.h> // size_t

typedef struct tridiagonal_solver_internal_t tridiagonal_solver_internal_t;

typedef struct {
  // size of linear system
  size_t nitems;
  // repeat several times
  size_t repeat_for;
  // for internal use, opaque pointer
  tridiagonal_solver_internal_t * internal;
} tridiagonal_solver_plan_t;

extern int tridiagonal_solver_init_plan(
    const size_t nitems,
    const size_t repeat_for,
    tridiagonal_solver_plan_t ** const tridiagonal_solver_plan
);

extern int tridiagonal_solver_exec(
    tridiagonal_solver_plan_t * const tridiagonal_solver_plan,
    // tri-diagonal matrix, lower, main, upper-diagonals
    const double * const l,
    const double * const m,
    const double * const u,
    // offset for main-diagonal components, can vary for each repeat
    const double * const m_offsets,
    // input and output
    double * const q
);

extern int tridiagonal_solver_destroy_plan(
    tridiagonal_solver_plan_t ** const tridiagonal_solver_plan
);

#endif // TRIDIAGONAL_SOLVER_H
