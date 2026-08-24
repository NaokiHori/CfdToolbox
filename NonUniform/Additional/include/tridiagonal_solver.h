#if !defined(TRIDIAGONAL_SOLVER_H)
#define TRIDIAGONAL_SOLVER_H

#include <stddef.h> // size_t

typedef struct {
  double lower;
  double main;
  double upper;
} tridiagonal_t;

// solve
//   (              prefactor * lower_diagonal) * p[i - 1]
// + (main_offset + prefactor *  main_diagonal) * p[i    ]
// + (              prefactor * upper_diagonal) * p[i + 1]
// = q[i]
extern int tridiagonal_solver_solve(
    const size_t nitems,
    // tridiagonal matrix
    // NOTE: the following elements are not used
    // - diagonals[0]->lower
    // - diagonals[nitems - 1]->upper
    const tridiagonal_t * const diagonals,
    // offset of main-diagonal component
    const double main_offset,
    // pre-factors in front of diagonals
    const double prefactor,
    // buffer used internally as a workspace
    double * const buffer,
    // input and output
    double * const q
);

#endif // TRIDIAGONAL_SOLVER_H
