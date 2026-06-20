#if !defined(RUNGE_KUTTA_H)
#define RUNGE_KUTTA_H

// runge-kutta integrator to solve odes:
// dp_i / dt + c_i * p_i = q_i(p_i)

#include <stddef.h> // size_t

typedef struct runge_kutta_t runge_kutta_t;

extern const size_t runge_kutta_step_max;

extern int initialize_runge_kutta(
    const size_t nitems,
    const double * const coefficients,
    runge_kutta_t ** const runge_kutta
);

extern int finalize_runge_kutta(
    runge_kutta_t ** const runge_kutta
);

extern int runge_kutta_set_time_step_size(
    runge_kutta_t * const runge_kutta,
    const double dt
);

extern int runge_kutta_update(
    runge_kutta_t * const runge_kutta,
    const double complex * const q,
    double complex * const p
);

#endif // RUNGE_KUTTA_H
