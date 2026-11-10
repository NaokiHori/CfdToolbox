#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "runge_kutta.h"

// Williamson, J. Comput. Phys. (35), 1980, along with integrating-factor technique

const size_t runge_kutta_step_max = 3;

struct runge_kutta_t {
  size_t step;
  double dt;
  bool is_dt_set;
  size_t nitems;
  const double * coefficients;
  double complex * gradients;
};

int initialize_runge_kutta(
    const size_t nitems,
    const double * const coefficients,
    runge_kutta_t ** const runge_kutta
) {
  *runge_kutta = malloc(sizeof(runge_kutta_t));
  (*runge_kutta)->step = 0;
  (*runge_kutta)->nitems = nitems;
  (*runge_kutta)->coefficients = coefficients;
  (*runge_kutta)->gradients = malloc(nitems * sizeof(double complex));
  return 0;
}

int finalize_runge_kutta(
    runge_kutta_t ** const runge_kutta
) {
  free(*runge_kutta);
  return 0;
}

int runge_kutta_set_time_step_size(
    runge_kutta_t * const runge_kutta,
    const double dt
) {
  if (0 != runge_kutta->step) {
    return 1;
  }
  runge_kutta->dt = dt;
  runge_kutta->is_dt_set = true;
  return 0;
}

int runge_kutta_update(
    runge_kutta_t * const runge_kutta,
    const double complex * const q,
    double complex * const p
) {
  bool * const is_dt_set = &runge_kutta->is_dt_set;
  if (!*is_dt_set) {
    return 1;
  }
  size_t * const step = &runge_kutta->step;
  const double dt = runge_kutta->dt;
  const size_t nitems = runge_kutta->nitems;
  const double * const c = runge_kutta->coefficients;
  double complex * const g = runge_kutta->gradients;
  if (0 == *step) {
#pragma omp parallel for
    for (size_t n = 0; n < nitems; n++) {
      g[n] = q[n] * dt;
      p[n] += (1. / 3.) * g[n];
    }
  } else if (1 == *step) {
#pragma omp parallel for
    for (size_t n = 0; n < nitems; n++) {
      g[n] = (- 5. / 9.) * g[n] + q[n] * dt;
      p[n] += (15. / 16.) * g[n];
    }
  } else {
#pragma omp parallel for
    for (size_t n = 0; n < nitems; n++) {
      g[n] = (- 153. / 128.) * g[n] + q[n] * dt;
      p[n] += (8. / 15.) * g[n];
    }
  }
  const double factor =
    0 == *step ? 1. / 3.
    : 1 == *step ? 5. / 12.
    : 1. / 4.;
#pragma omp parallel for
  for (size_t n = 0; n < nitems; n++) {
    const double exponential = exp(- c[n] * factor * dt);
    p[n] *= exponential;
    g[n] *= exponential;
  }
  *step = (*step + 1) % runge_kutta_step_max;
  if (0 == *step) {
    // time-step size should be set before the next integration starts
    *is_dt_set = false;
  }
  return 0;
}

#if defined(TEST_RUNGE_KUTTA)

#include <stdio.h>

// solve dp / dt + c * p = q(p)

static double compute_right_hand_side(
    const double p
) {
  return - p * p;
}

static double complex compute_exact_solution(
    const double coefficient,
    const double time
) {
  const double exponential = exp(- coefficient * time);
  return coefficient * exponential / (coefficient + 1. - exponential);
}

int main(
    void
) {
  const double time_max = 1.;
#define NITEMS 5
  for (size_t step_max = (1 << 3); step_max < (1 << 12); step_max <<= 1) {
    const double dt = time_max / step_max;
    const double coefficients[NITEMS] = {1e-3, 1e-2, 1e-1, 1e+0, 1e+1};
    double complex p[NITEMS] = {1., 1., 1., 1., 1.};
    double complex q[NITEMS] = {0.};
    runge_kutta_t * runge_kutta = NULL;
    initialize_runge_kutta(NITEMS, coefficients, &runge_kutta);
    for (size_t step = 0; step < step_max; step++) {
      runge_kutta_set_time_step_size(runge_kutta, dt);
      // step 0
      for (size_t n = 0; n < NITEMS; n++) {
        q[n] = compute_right_hand_side(p[n]);
      }
      runge_kutta_update(runge_kutta, q, p);
      // step 1
      for (size_t n = 0; n < NITEMS; n++) {
        q[n] = compute_right_hand_side(p[n]);
      }
      runge_kutta_update(runge_kutta, q, p);
      // step 2
      for (size_t n = 0; n < NITEMS; n++) {
        q[n] = compute_right_hand_side(p[n]);
      }
      runge_kutta_update(runge_kutta, q, p);
    }
    finalize_runge_kutta(&runge_kutta);
    for (size_t n = 0; n < NITEMS; n++) {
      printf(
          "%zu % .7e % .7e\n",
          n,
          dt,
          cabs(p[n] - compute_exact_solution(coefficients[n], time_max))
      );
    }
  }
  return 0;
}

#endif // TEST_RUNGE_KUTTA

