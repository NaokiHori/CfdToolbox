#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stddef.h> // size_t

typedef struct {
  double lx;
  double ly;
  size_t nx;
  size_t ny;
  double dx;
  double dy;
} domain_t;

extern int domain_init(
    domain_t * const domain
);

#endif // DOMAIN_H
