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
    const double * const l,
    const double * const m,
    const double * const u,
    const double m_offset,
    double * const v,
    double * const q
) {
  // divide the first row by center-diagonal term
  v[0] = u[0] / (m[0] + m_offset);
  q[0] = q[0] / (m[0] + m_offset);
  // forward sweep
  for (size_t i = 1; i < nitems - 1; i++) {
    // assume positive-definite system
    //   to skip zero-division checks
    const double val = 1. / (m[i] + m_offset - l[i] * v[i - 1]);
    v[i] = val * u[i];
    q[i] = val * (q[i] - l[i] * q[i - 1]);
  }
  // last row, do the same thing but consider singularity (degeneracy)
  const double val = m[nitems - 1] + m_offset - l[nitems - 1] * v[nitems - 2];
  if (DBL_EPSILON < myfabs(val)) {
    q[nitems - 1] = 1. / val * (q[nitems - 1] - l[nitems - 1] * q[nitems - 2]);
  } else {
    // singular
    q[nitems - 1] = 0.;
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

