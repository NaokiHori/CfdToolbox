#if !defined(TIME_MARCHER_H)
#define TIME_MARCHER_H

// specify implicit time-marching scheme
// - 0.0: Euler-explicit
// - 0.5: Crank-Nicolson
// - 1.0: Euler-implicit
// other values may lead to undefined behavior
#define IMPLICIT_FACTOR 0.5

typedef struct {
  double ** non_linear;
  double ** increment;
} runge_kutta_buffer_t;

// use three-stage Runge-Kutta
#define RUNGE_KUTTA_STEP_MAX 3

typedef struct {
  double alpha[RUNGE_KUTTA_STEP_MAX];
  double beta[RUNGE_KUTTA_STEP_MAX];
  double gamma[RUNGE_KUTTA_STEP_MAX];
} runge_kutta_coefs_t;

extern const runge_kutta_coefs_t RUNGE_KUTTA_COEFS;

#endif // TIME_MARCHER_H
