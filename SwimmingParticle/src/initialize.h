#if !defined(INITIALIZE_H)
#define INITIALIZE_H

#include "domain.h"

extern int initialize(
    const domain_t * const domain,
    double ** const concentration
);

#endif // INITIALIZE_H
