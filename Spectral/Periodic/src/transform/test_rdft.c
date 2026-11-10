#if defined(TEST_RDFT)

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rdft.h"

static const double pi = 3.141592653589793238462643383;

static int rdft_naive_exec_f(
    const size_t nitems,
    double * const xs
) {
  double complex * const ys = malloc(nitems * sizeof(double complex));
  double complex * const zs = malloc(nitems * sizeof(double complex));
  for (size_t n = 0; n < nitems; n++) {
    ys[n] = xs[n] + 0. * I;
  }
  for (size_t k = 0; k < nitems; k++) {
    zs[k] = 0.;
    for (size_t n = 0; n < nitems; n++) {
      zs[k] += ys[n] * cexp(- 2. * pi * n * k * I / nitems);
    }
  }
  for (size_t k = 0; k < nitems; k++) {
    xs[k] = k < nitems / 2 + 1 ? creal(zs[k]) : cimag(zs[k]);
  }
  free(ys);
  free(zs);
  return 0;
}

static int rdft_naive_exec_b(
    const size_t nitems,
    double * const xs
) {
  double complex * const ys = malloc(nitems * sizeof(double complex));
  double complex * const zs = malloc(nitems * sizeof(double complex));
  for (size_t k = 0; k < nitems; k++) {
    if (k == 0) {
      ys[k] = xs[k];
    } else if (k < nitems / 2) {
      ys[k] = xs[k] - xs[nitems - k] * I;
    } else if (k == nitems / 2) {
      ys[k] = xs[k];
    } else {
      ys[k] = xs[nitems - k] + xs[k] * I;
    }
  }
  for (size_t n = 0; n < nitems; n++) {
    zs[n] = 0.;
    for (size_t k = 0; k < nitems; k++) {
      zs[n] += ys[k] * cexp(+ 2. * pi * n * k * I / nitems);
    }
  }
  for (size_t n = 0; n < nitems; n++) {
    xs[n] = creal(zs[n]) / nitems;
  }
  free(ys);
  free(zs);
  return 0;
}

static int set_values(
    const size_t nitems,
    double * const xs,
    double * const ys
) {
  for (size_t i = 0; i < nitems; i++) {
    xs[i] = sin(2. * pi * i / nitems) + 1. * rand() / RAND_MAX;
    ys[i] = xs[i];
  }
  return 0;
}

static int test_0_in_place(
    const size_t nitems,
    rdft_plan_t * const plan
) {
  double * const xs = malloc(nitems * sizeof(double));
  double * const ys = malloc(nitems * sizeof(double));
  set_values(nitems, xs, ys);
  rdft_exec_f(plan, xs, xs);
  rdft_exec_b(plan, xs, xs);
  double error = 0.;
  for (size_t i = 0; i < nitems; i++) {
    error += pow(xs[i] - ys[i], 2.);
  }
  error = sqrt(error / nitems);
  printf("residual: % .7e\n", error);
  free(xs);
  free(ys);
  return 0;
}

static int test_0_out_of_place(
    const size_t nitems,
    rdft_plan_t * const plan
) {
  double * const xs = malloc(nitems * sizeof(double));
  double * const ys = malloc(nitems * sizeof(double));
  double * const zs = malloc(nitems * sizeof(double));
  set_values(nitems, xs, ys);
  rdft_exec_f(plan, xs, zs);
  rdft_exec_b(plan, zs, xs);
  double error = 0.;
  for (size_t i = 0; i < nitems; i++) {
    error += pow(xs[i] - ys[i], 2.);
  }
  error = sqrt(error / nitems);
  printf("residual: % .7e\n", error);
  free(xs);
  free(ys);
  free(zs);
  return 0;
}

static int test_1_in_place(
    const size_t nitems,
    rdft_plan_t * const plan
) {
  double * const xs = malloc(nitems * sizeof(double));
  double * const ys = malloc(nitems * sizeof(double));
  set_values(nitems, xs, ys);
  rdft_exec_f(plan, xs, xs);
  rdft_naive_exec_f(nitems, ys);
  double error = 0.;
  for (size_t i = 0; i < nitems; i++) {
    error += pow(xs[i] - ys[i], 2.);
  }
  error = sqrt(error / nitems);
  printf("residual: % .7e\n", error);
  free(xs);
  free(ys);
  return 0;
}

static int test_1_out_of_place(
    const size_t nitems,
    rdft_plan_t * const plan
) {
  double * const xs = malloc(nitems * sizeof(double));
  double * const ys = malloc(nitems * sizeof(double));
  double * const zs = malloc(nitems * sizeof(double));
  set_values(nitems, xs, ys);
  rdft_exec_f(plan, xs, zs);
  rdft_naive_exec_f(nitems, ys);
  double error = 0.;
  for (size_t i = 0; i < nitems; i++) {
    error += pow(zs[i] - ys[i], 2.);
  }
  error = sqrt(error / nitems);
  printf("residual: % .7e\n", error);
  free(xs);
  free(ys);
  free(zs);
  return 0;
}

static int test_2_in_place(
    const size_t nitems,
    rdft_plan_t * const plan
) {
  double * const xs = malloc(nitems * sizeof(double));
  double * const ys = malloc(nitems * sizeof(double));
  set_values(nitems, xs, ys);
  rdft_exec_b(plan, xs, xs);
  rdft_naive_exec_b(nitems, ys);
  double error = 0.;
  for (size_t i = 0; i < nitems; i++) {
    error += pow(xs[i] - ys[i], 2.);
  }
  error = sqrt(error / nitems);
  printf("residual: % .7e\n", error);
  free(xs);
  free(ys);
  return 0;
}

static int test_2_out_of_place(
    const size_t nitems,
    rdft_plan_t * const plan
) {
  double * const xs = malloc(nitems * sizeof(double));
  double * const ys = malloc(nitems * sizeof(double));
  double * const zs = malloc(nitems * sizeof(double));
  set_values(nitems, xs, ys);
  rdft_exec_b(plan, xs, zs);
  rdft_naive_exec_b(nitems, ys);
  double error = 0.;
  for (size_t i = 0; i < nitems; i++) {
    error += pow(zs[i] - ys[i], 2.);
  }
  error = sqrt(error / nitems);
  printf("residual: % .7e\n", error);
  free(xs);
  free(ys);
  free(zs);
  return 0;
}

int main(
    void
) {
  const size_t nitems = 512;
  rdft_plan_t * plan = NULL;
  if (0 != rdft_init_plan(nitems, &plan)) {
    return 1;
  }
  test_0_in_place(nitems, plan);
  test_0_out_of_place(nitems, plan);
  test_1_in_place(nitems, plan);
  test_1_out_of_place(nitems, plan);
  test_2_in_place(nitems, plan);
  test_2_out_of_place(nitems, plan);
  if (0 != rdft_destroy_plan(&plan)) {
    return 1;
  }
  return 0;
}

#endif // TEST_RDFT

extern char dummy;

