#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stddef.h>

typedef struct {
  double lx;
  double ly;
  size_t nx;
  size_t ny;
  size_t knx;
  size_t kny;
  double * kxs;
  double * kys;
} domain_t;

extern int initialize_domain(
    const double lx,
    const double ly,
    const size_t nx,
    const size_t ny,
    domain_t * const domain
);

extern int finalize_domain(
    domain_t * const domain
);

#endif // DOMAIN_H
