#if !defined(TRIDIAGONAL_SOLVER_H)
#define TRIDIAGONAL_SOLVER_H

#include <stddef.h> // size_t

extern int tridiagonal_solver_solve(
    const size_t nitems,
    // tridiagonal matrix
    const double * const lower_diagonals,
    const double * const main_diagonals,
    const double * const upper_diagonals,
    // offset of main-diagonal component
    const double main_offset,
    // buffer used internally as a workspace
    double * const buffer,
    // input and output
    double * const q
);

#endif // TRIDIAGONAL_SOLVER_H
