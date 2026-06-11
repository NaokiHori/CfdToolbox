#if !defined(DECIDE_TIME_STEP_SIZE_H)
#define DECIDE_TIME_STEP_SIZE_H

#include "domain.h"

extern int decide_time_step_size(
    const double peclet,
    const domain_t * const domain,
    double ** const stream_function,
    double * const dt
);

#endif // DECIDE_TIME_STEP_SIZE_H
