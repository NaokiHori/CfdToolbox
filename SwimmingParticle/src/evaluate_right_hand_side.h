#if !defined(EVALUATE_RIGHT_HAND_SIDE_H)
#define EVALUATE_RIGHT_HAND_SIDE_H

#include <stddef.h>
#include "domain.h"
#include "time_marcher.h"

extern int evaluate_right_hand_side(
    const double peclet,
    const domain_t * const domain,
    const size_t runge_kutta_step,
    double ** const concentration,
    double ** const stream_function,
    runge_kutta_buffer_t * const runge_kutta_buffer
);

#endif // EVALUATE_RIGHT_HAND_SIDE_H
