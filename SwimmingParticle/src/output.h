#if !defined(OUTPUT_H)
#define OUTPUT_H

#include "domain.h"

extern int output(
    const domain_t * const domain,
    double ** const concentration
);

#endif // OUTPUT_H
