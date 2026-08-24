#include <stdio.h>
#include <float.h> // DBL_EPSILON
#include "memory.h"
#include "tridiagonal_solver.h"

static double myfabs(
    const double v
) {
  return v < 0. ? - v : v;
}

int tridiagonal_solver_solve(
    const size_t nitems,
    const tridiagonal_t * const diagonals,
    const double prefactor,
    const double m_offset,
    double * const v,
    double * const q
) {
  // divide the first row by center-diagonal term
  {
    const tridiagonal_t * const diagonal = diagonals + 0;
    v[0] = prefactor * diagonal->upper / (prefactor * diagonal->main + m_offset);
    q[0] = q[0] / (prefactor * diagonal->main + m_offset);
  }
  // forward sweep
  for (size_t i = 1; i < nitems - 1; i++) {
    const tridiagonal_t * const diagonal = diagonals + i;
    // assume positive-definite system
    //   to skip zero-division checks
    const double value = 1. / (
        + prefactor * diagonal->main
        + m_offset
        - prefactor * diagonal->lower * v[i - 1]
    );
    v[i] = value * prefactor * diagonal->upper;
    q[i] = value * (q[i] - prefactor * diagonal->lower * q[i - 1]);
  }
  // last row, do the same thing but consider singularity (degeneracy)
  {
    const tridiagonal_t * const diagonal = diagonals + nitems - 1;
    const double value =
      + prefactor * diagonal->main
      + m_offset
      - prefactor * diagonal->lower * v[nitems - 2];
    if (DBL_EPSILON < myfabs(value)) {
      q[nitems - 1] = 1. / value * (q[nitems - 1] - prefactor * diagonal->lower * q[nitems - 2]);
    } else {
      // singular
      q[nitems - 1] = 0.;
    }
  }
  // backward substitution
  for (size_t i = nitems - 2; ; i--) {
    q[i] -= v[i] * q[i + 1];
    if (0 == i) {
      break;
    }
  }
  return 0;
}

