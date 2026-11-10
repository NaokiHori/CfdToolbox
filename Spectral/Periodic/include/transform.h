#if !defined(TRANSFORM_H)
#define TRANSFORM_H

#include <complex.h>
#include <stddef.h>
#include "domain.h"

typedef struct transformer_t transformer_t;

extern int initialize_transformer(
    const size_t nx,
    const size_t ny,
    transformer_t ** const transformer
);

extern int finalize_transformer(
    transformer_t ** const transformer
);

extern int transform_physical_to_spectral(
    transformer_t * const transformer,
    double * const physical,
    double complex * const spectral
);

extern int transform_spectral_to_physical(
    transformer_t * const transformer,
    const double complex * const spectral,
    double * const physical
);

#endif // TRANSFORM_H
