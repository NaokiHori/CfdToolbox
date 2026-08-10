#include <stdio.h>
#include <math.h> // fmin, fmax, fabs, pow
#include "logger.h"
#include "param.h"
#include "./decide_dt.h"

static const struct {
  double adv;
  double dif;
} safety_factors = {
  .adv = 0.25,
  .dif = 0.95,
};

static int decide_dt_adv(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double small = 1.e-8;
  double ** const ux = flow_field->ux;
  double ** const uy = flow_field->uy;
  *dt = 1.;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      const double denominator = fmax(small, fabs(ux[j][i]));
      *dt = fmin(*dt, x_face_scale_factors[i] / denominator);
    }
  }
  for (size_t j = 2; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double denominator = fmax(small, fabs(uy[j][i]));
      *dt = fmin(*dt, y_face_scale_factors[j] / denominator);
    }
  }
  *dt *= safety_factors.adv;
  return 0;
}

static int decide_dt_dif(
    const domain_t * const domain,
    double * const dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const coordinate_t * const y_coordinate = &domain->y_coordinate;
  const double * const x_center_scale_factors = x_coordinate->center_scale_factors;
  const double * const y_center_scale_factors = y_coordinate->center_scale_factors;
  double dx = x_coordinate->length;
  double dy = y_coordinate->length;
  for (size_t i = 1; i <= nx; i++) {
    dx = fmin(dx, x_center_scale_factors[i]);
  }
  for (size_t j = 1; j <= ny; j++) {
    dy = fmin(dy, y_center_scale_factors[j]);
  }
  *dt = Re * 0.5 / NDIMS * pow(fmin(dx, dy), 2.);
  *dt *= safety_factors.dif;
  return 0;
}

int decide_dt(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  double dt_adv = 0.;
  double dt_dif = 0.;
  if (0 != decide_dt_adv(domain, flow_field, &dt_adv)) {
    LOGGER_FAILURE("failed to find advective time-step constraint");
    goto abort;
  }
  if (0 != decide_dt_dif(domain, &dt_dif)) {
    LOGGER_FAILURE("failed to find diffusive time-step constraint");
    goto abort;
  }
  *dt = fmin(dt_adv, dt_dif);
  return 0;
abort:
  LOGGER_FAILURE("failed to find time-step size");
  return 1;
}

