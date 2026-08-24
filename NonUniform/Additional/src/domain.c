#include <math.h>
#include "domain.h"
#include "memory.h"
#include "param.h"

// grid-streching indicator
// 0: extremely stretched
// 1: uniform
static int build_coordinate(
    const double length,
    const size_t nitems,
    const double grad,
    coordinate_t * const coordinate
) {
  coordinate->length = length;
  // coordinates
  double ** const faces = &coordinate->faces;
  double ** const centers = &coordinate->centers;
  *faces = memory_alloc(nitems + 2, sizeof(double));
  *centers = memory_alloc(nitems + 2, sizeof(double));
  (*faces)[0] = nan("");
  (*faces)[1] = 0.;
  for (size_t i = 2; i < nitems + 1; i++) {
    double * const face = *faces + i;
    const double x = 1. * (i - 1) / nitems;
    *face =
      + (2. * grad - 2.) * (x * x * x)
      + (- 3. * grad + 3.) * (x * x)
      + grad * x;
    *face *= length;
  }
  (*faces)[nitems + 1] = length;
  (*centers)[0] = (*faces)[1];
  for (size_t i = 1; i < nitems + 1; i++) {
    (*centers)[i] = 0.5 * (*faces)[i] + 0.5 * (*faces)[i + 1];
  }
  (*centers)[nitems + 1] = (*faces)[nitems + 1];
  // scale factors
  double ** const face_scale_factors = &coordinate->face_scale_factors;
  double ** const center_scale_factors = &coordinate->center_scale_factors;
  *face_scale_factors = memory_alloc(nitems + 2, sizeof(double));
  *center_scale_factors = memory_alloc(nitems + 2, sizeof(double));
  (*face_scale_factors)[0] = nan("");
  for (size_t i = 1; i < nitems + 2; i++) {
    (*face_scale_factors)[i] = (*centers)[i] - (*centers)[i - 1];
  }
  (*center_scale_factors)[0] = nan("");
  for (size_t i = 1; i < nitems + 1; i++) {
    (*center_scale_factors)[i] = (*faces)[i + 1] - (*faces)[i];
  }
  (*center_scale_factors)[nitems + 1] = nan("");
  // laplace operators
  const tridiagonal_t operator_at_not_in_use = {
    .lower = nan(""),
    .main = nan(""),
    .upper = nan(""),
  };
  tridiagonal_t ** const face_laplace_operators = &coordinate->face_laplace_operators;
  *face_laplace_operators = memory_alloc(nitems + 2, sizeof(tridiagonal_t));
  (*face_laplace_operators)[0] = operator_at_not_in_use;
  (*face_laplace_operators)[1] = operator_at_not_in_use;
  for (size_t i = 2; i < nitems + 1; i++) {
    tridiagonal_t * const face_laplace_operator = *face_laplace_operators + i;
    double * const lower = &face_laplace_operator->lower;
    double * const main = &face_laplace_operator->main;
    double * const upper = &face_laplace_operator->upper;
    *lower = 1. / (*face_scale_factors)[i] / (*center_scale_factors)[i - 1];
    *upper = 1. / (*face_scale_factors)[i] / (*center_scale_factors)[i    ];
    *main = - *lower - *upper;
  }
  (*face_laplace_operators)[nitems + 1] = operator_at_not_in_use;
  tridiagonal_t ** const center_laplace_operators = &coordinate->center_laplace_operators;
  *center_laplace_operators = memory_alloc(nitems + 2, sizeof(tridiagonal_t));
  (*center_laplace_operators)[0] = operator_at_not_in_use;
  for (size_t i = 1; i < nitems + 1; i++) {
    tridiagonal_t * const center_laplace_operator = *center_laplace_operators + i;
    double * const lower = &center_laplace_operator->lower;
    double * const main = &center_laplace_operator->main;
    double * const upper = &center_laplace_operator->upper;
    *lower = 1. / (*center_scale_factors)[i] / (*face_scale_factors)[i    ];
    *upper = 1. / (*center_scale_factors)[i] / (*face_scale_factors)[i + 1];
    *main = - *lower - *upper;
  }
  (*center_laplace_operators)[nitems + 1] = operator_at_not_in_use;
  return 0;
}

static int coordinate_finalize(
    coordinate_t * const coordinate
) {
  memory_free(coordinate->faces);
  memory_free(coordinate->centers);
  memory_free(coordinate->face_scale_factors);
  memory_free(coordinate->center_scale_factors);
  memory_free(coordinate->face_laplace_operators);
  memory_free(coordinate->center_laplace_operators);
  return 0;
}

int domain_init(
    domain_t * const domain
) {
  const double lengths[NDIMS] = {0.25, 1.};
  const size_t nitems[NDIMS] = {NX, NY};
  if (0 != build_coordinate(
      lengths[0],
      nitems[0],
      0.25,
      &domain->x_coordinate
  )) {
    return 1;
  }
  if (0 != build_coordinate(
      lengths[1],
      nitems[1],
      0.125,
      &domain->y_coordinate
  )) {
    return 1;
  }
  return 0;
}

int domain_finalize(
    domain_t * const domain
) {
  coordinate_finalize(&domain->x_coordinate);
  coordinate_finalize(&domain->y_coordinate);
  return 0;
}

