#if !defined(DOMAIN_H)
#define DOMAIN_H

#include "tridiagonal_solver.h"

#define NDIMS 2

#define NX 33
#define NY 129

typedef struct {
  double length;
  double * faces;
  double * centers;
  double * face_scale_factors;
  double * center_scale_factors;
  tridiagonal_t * face_laplace_operators;
  tridiagonal_t * center_laplace_operators;
} coordinate_t;

typedef struct {
  coordinate_t x_coordinate;
  coordinate_t y_coordinate;
} domain_t;

extern int domain_init(
    domain_t * const domain
);

extern int domain_finalize(
    domain_t * const domain
);

#endif // DOMAIN_H
