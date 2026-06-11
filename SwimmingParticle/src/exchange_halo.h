#if !defined(EXCHANGE_HALO_H)
#define EXCHANGE_HALO_H

#include "domain.h"

extern int exchange_halo(
    const domain_t * const domain,
    double ** const concentration
);

#endif // EXCHANGE_HALO_H
