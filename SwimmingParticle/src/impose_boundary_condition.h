#if !defined(IMPOSE_BOUNDARY_CONDITION_H)
#define IMPOSE_BOUNDARY_CONDITION_H

#include "domain.h"

extern int impose_boundary_condition(
    const domain_t * const domain,
    double ** const concentration
);

#endif // IMPOSE_BOUNDARY_CONDITION_H
