#if defined(TRIDIAGONAL_SOLVER_TEST)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory.h"
#include "tridiagonal_solver.h"

#define MY_ASSERT(cond) \
  if (!(cond)) { \
    fprintf(stderr, "Test failed: %s (%d)\n", __func__, __LINE__); \
    retval += 1; \
  }

static const double small = 1.e-14;

int main (
    void
) {
  int retval = 0;
  // solve
  // -2  1  0  0 | -3
  //  1 -2  1  0 | -1
  //  0  1 -2  1 | +1
  //  0  0  1 -2 | +3
  // answer: (+2, +1, -1, -2)
  const size_t nitems = 4;
  double * l = memory_alloc(nitems, sizeof(double));
  double * m = memory_alloc(nitems, sizeof(double));
  double * u = memory_alloc(nitems, sizeof(double));
  double * q = memory_alloc(nitems, sizeof(double));
  double * buffer = memory_alloc(nitems, sizeof(double));
  for (size_t i = 0; i < nitems; i++) {
    l[i] = + 1.;
    m[i] = - 2.;
    u[i] = + 1.;
  }
  q[0] = - 3.;
  q[1] = - 1.;
  q[2] = + 1.;
  q[3] = + 3.;
  MY_ASSERT(0 == tridiagonal_solver_solve(nitems, l, m, u, 0., buffer, q));
  MY_ASSERT(fabs(q[0] - 2.) < small);
  MY_ASSERT(fabs(q[1] - 1.) < small);
  MY_ASSERT(fabs(q[2] + 1.) < small);
  MY_ASSERT(fabs(q[3] + 2.) < small);
  memory_free(l);
  memory_free(m);
  memory_free(u);
  memory_free(q);
  memory_free(buffer);
  return retval;
}

#else
extern char dummy;
#endif // TRIDIAGONAL_SOLVER_TEST
