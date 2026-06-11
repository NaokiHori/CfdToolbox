#if !defined(FLUID_H)
#define FLUID_H

#include <stddef.h>
#include "domain.h"

extern int initialize_fluid_distribution_function(
    df_t * const f
);

extern int compute_fluid_macroscopic_field(
    const df_t * const f,
    scalar_t * const density,
    vector_t * const velocity
);

extern int process_fluid_collision(
    const double diffusivity,
    const double acceleration,
    const scalar_t * const density,
    const vector_t * const velocity,
    const scalar_t * const temperature,
    df_t * const f
);

extern int process_fluid_streaming(
    df_t * const f
);

#endif // FLUID_H
