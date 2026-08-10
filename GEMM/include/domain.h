#if !defined(DOMAIN_H)
#define DOMAIN_H

#include <stddef.h> // size_t

typedef struct {
  double length;
  double * faces;
  double * centers;
  double * face_scale_factors;
  double * center_scale_factors;
} coordinate_t;

typedef struct {
  coordinate_t x_coordinate;
  coordinate_t y_coordinate;
  size_t nx;
  size_t ny;
} domain_t;

extern int domain_init(
    domain_t * const domain
);

extern int domain_finalize(
    domain_t * const domain
);

#endif // DOMAIN_H
