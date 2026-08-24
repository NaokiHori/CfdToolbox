#if !defined(RUNGE_KUTTA_H)
#define RUNGE_KUTTA_H

#include <stddef.h>

extern const size_t N_RUNGE_KUTTA_STEP;

typedef struct {
  size_t step;
  double explicit_previous;
  double explicit_current;
  double implicit;
} runge_kutta_coefficients_t;

extern int runge_kutta_coefficients_setup(
    const size_t step,
    runge_kutta_coefficients_t * const coefficients
);

#endif // RUNGE_KUTTA_H
