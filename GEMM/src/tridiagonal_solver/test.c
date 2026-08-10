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

#define REPORT_SUCCESS(objective) \
  if (0 == retval) { \
    fprintf(stderr, "Test passed: %s (%s)\n", __func__, objective); \
  }

#define REPEAT_FOR 2

static const double small = 1.e-14;

static int test0 (
    void
) {
  int retval = 0;
  // solve
  // -2  1  0  0 | (-3, -6)
  //  1 -2  1  0 | (-1, -2)
  //  0  1 -2  1 | (+1, +2)
  //  0  0  1 -2 | (+3, +6)
  // answer: (+2, +1, -1, -2), (+4, +2, -2, -4)
  const size_t nitems = 4;
  const char objective[] = "non-singular";
  double * l = memory_alloc(nitems, sizeof(double));
  double * c = memory_alloc(nitems, sizeof(double));
  double * u = memory_alloc(nitems, sizeof(double));
  double * q = memory_alloc(nitems * REPEAT_FOR, sizeof(double));
  double * x = memory_alloc(nitems * REPEAT_FOR, sizeof(double));
  for (size_t i = 0; i < nitems; i++) {
    l[i] = + 1.;
    c[i] = - 2.;
    u[i] = + 1.;
  }
  q[0] = - 3.;
  q[1] = - 1.;
  q[2] = + 1.;
  q[3] = + 3.;
  q[4] = - 6.;
  q[5] = - 2.;
  q[6] = + 2.;
  q[7] = + 6.;
  for (size_t i = 0; i < nitems * REPEAT_FOR; i++) {
    x[i] = q[i];
  }
  tridiagonal_solver_plan_t * plan = NULL;
  MY_ASSERT(0 == tridiagonal_solver_init_plan(nitems, REPEAT_FOR, &plan));
  MY_ASSERT(NULL != plan);
  MY_ASSERT(0 == tridiagonal_solver_exec(plan, l, c, u, (double [REPEAT_FOR]){0., 0.}, x));
  MY_ASSERT(fabs(x[0] - 2.) < small);
  MY_ASSERT(fabs(x[1] - 1.) < small);
  MY_ASSERT(fabs(x[2] + 1.) < small);
  MY_ASSERT(fabs(x[3] + 2.) < small);
  MY_ASSERT(fabs(x[4] - 4.) < small);
  MY_ASSERT(fabs(x[5] - 2.) < small);
  MY_ASSERT(fabs(x[6] + 2.) < small);
  MY_ASSERT(fabs(x[7] + 4.) < small);
  MY_ASSERT(fabs(            + c[0] * x[0] + u[0] * x[1] - q[0]) < small);
  MY_ASSERT(fabs(l[1] * x[0] + c[1] * x[1] + u[1] * x[2] - q[1]) < small);
  MY_ASSERT(fabs(l[2] * x[1] + c[2] * x[2] + u[2] * x[3] - q[2]) < small);
  MY_ASSERT(fabs(l[3] * x[2] + c[3] * x[3]               - q[3]) < small);
  MY_ASSERT(fabs(            + c[0] * x[4] + u[0] * x[5] - q[4]) < small);
  MY_ASSERT(fabs(l[1] * x[4] + c[1] * x[5] + u[1] * x[6] - q[5]) < small);
  MY_ASSERT(fabs(l[2] * x[5] + c[2] * x[6] + u[2] * x[7] - q[6]) < small);
  MY_ASSERT(fabs(l[3] * x[6] + c[3] * x[7]               - q[7]) < small);
  MY_ASSERT(0 == tridiagonal_solver_destroy_plan(&plan));
  MY_ASSERT(NULL == plan);
  memory_free(l);
  memory_free(c);
  memory_free(u);
  memory_free(q);
  memory_free(x);
  REPORT_SUCCESS(objective);
  return retval;
}

int main (
    void
) {
  int retval = 0;
  retval += test0();
  return retval;
}

#else
extern char dummy;
#endif // TRIDIAGONAL_SOLVER_TEST
