#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stddef.h>

typedef struct {
  // number of grids
  size_t nx;
  size_t ny;
  // radial cell face / center positions
  double * xf;
  double * xc;
  // radial scale factors at radial cell face / center
  double * hxxf;
  double * hxxc;
  // azimuthal scale factors at radial cell face / center
  double * hyxf;
  double * hyxc;
  // jacobian determinants at radial cell face / center
  double * jdxf;
  double * jdxc;
  // azimuthal grid size
  double dy;
} domain_t;

extern int domain_initialize(
    const size_t nx,
    const size_t ny,
    const double outer_radius,
    domain_t * const domain
);

extern int domain_finalize(
    domain_t * const domain
);

#endif // DOMAIN_H
