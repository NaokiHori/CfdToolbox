#if !defined(TEMPERATURE_H)
#define TEMPERATURE_H

#include <stddef.h>
#include "domain.h"

extern int initialize_temperature_distribution_function(
    df_t * const g
);

extern int compute_temperature_macroscopic_field(
    const df_t * const g,
    scalar_t * const temperature
);

extern int process_temperature_collision(
    const double diffusivity,
    const vector_t * const velocity,
    const scalar_t * const temperature,
    df_t * const g
);

extern int process_temperature_streaming(
    df_t * const g
);

#endif // TEMPERATURE_H
