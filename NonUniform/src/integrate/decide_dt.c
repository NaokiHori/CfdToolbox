#include <stdio.h>
#include <math.h> // fmin, fmax, fabs, pow
#include "logger.h"
#include "./decide_dt.h"

static const double dt_max = 1e-2;

static const struct {
  double advection;
  double buoyancy;
} safety_factors = {
  .advection = 0.25,
  .buoyancy = 0.25,
};

static int decide_dt_advection(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  const double * const x_face_scale_factors = domain->x_coordinate.face_scale_factors;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double small = 1.e-8;
  const array_t * const ux_array = flow_field->ux;
  const array_t * const uy_array = flow_field->uy;
  const double (* const ux)[NX + 2] = ux_array->buffer;
  const double (* const uy)[NX + 2] = uy_array->buffer;
  *dt = 1.;
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 2; i <= NX; i++) {
      const double denominator = fmax(small, fabs(ux[j][i]));
      *dt = fmin(*dt, x_face_scale_factors[i] / denominator);
    }
  }
  for (size_t j = 2; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double denominator = fmax(small, fabs(uy[j][i]));
      *dt = fmin(*dt, y_face_scale_factors[j] / denominator);
    }
  }
  *dt *= safety_factors.advection;
  return 0;
}

static int decide_dt_buoyancy(
    const parameters_t * const parameters,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  const double density_ratio = parameters->density_ratio;
  const double * const y_face_scale_factors = domain->y_coordinate.face_scale_factors;
  const double small = 1.e-8;
  const array_t * const temperature_array = flow_field->temperature;
  const array_t * const salinity_array = flow_field->salinity;
  const double (* const temperature)[NX + 2] = temperature_array->buffer;
  const double (* const salinity)[NX + 2] = salinity_array->buffer;
  *dt = 1.;
  for (size_t j = 2; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double acceleration =
        (0.5 * temperature[j - 1][i] + 0.5 * temperature[j][i])
        -
        density_ratio * (0.5 * salinity[j - 1][i] + 0.5 * salinity[j][i]);
      const double denominator = fmax(small, fabs(acceleration));
      *dt = fmin(*dt, sqrt(y_face_scale_factors[j] / denominator));
    }
  }
  *dt *= safety_factors.buoyancy;
  return 0;
}

int decide_dt(
    const parameters_t * const parameters,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const dt
) {
  double dt_advection = 0.;
  if (0 != decide_dt_advection(domain, flow_field, &dt_advection)) {
    goto abort;
  }
  double dt_buoyancy = 0.;
  if (0 != decide_dt_buoyancy(parameters, domain, flow_field, &dt_buoyancy)) {
    goto abort;
  }
  //
  *dt = dt_max;
  *dt = fmin(*dt, dt_advection);
  *dt = fmin(*dt, dt_buoyancy);
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

