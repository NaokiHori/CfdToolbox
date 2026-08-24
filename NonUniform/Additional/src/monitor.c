#include <stdio.h>
#include <errno.h>
#include <math.h>
#include "logger.h"
#include "flow_field.h"
#include "./monitor.h"

#define ROOT_DIRECTORY "output/log/"

static int output(
    const size_t step,
    const double time,
    const char file_name[],
    const size_t nitems,
    const double * quantities
) {
  errno = 0;
  FILE * const fp = fopen(file_name, "a");
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  fprintf(fp, "%10zu % .15e ", step, time);
  for (size_t n = 0; n < nitems; n++) {
    fprintf(fp, "% .15e%c", quantities[n], nitems - 1 == n ? '\n' : ' ');
  }
  fclose(fp);
  return 0;
}

static int monitor_divergence(
    const size_t step,
    const double time,
    const domain_t * const domain,
    const flow_field_t * const flow_field
) {
  const char file_name[] = ROOT_DIRECTORY "divergence.dat";
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const coordinate_t * const y_coordinate = &domain->y_coordinate;
  const double * const x_center_scale_factors = x_coordinate->center_scale_factors;
  const double * const y_center_scale_factors = y_coordinate->center_scale_factors;
  const array_t * const ux_array = flow_field->ux;
  const array_t * const uy_array = flow_field->uy;
  const double (* const ux)[NX + 2] = ux_array->buffer;
  const double (* const uy)[NX + 2] = uy_array->buffer;
  double div_max = 0.;
  double div_sum = 0.;
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double dux = - ux[j    ][i    ] + ux[j    ][i + 1];
      const double duy = - uy[j    ][i    ] + uy[j + 1][i    ];
      const double div =
        + 1. / x_center_scale_factors[i] * dux
        + 1. / y_center_scale_factors[j] * duy;
      div_max = fmax(div_max, fabs(div));
      div_sum = div_sum + div;
    }
  }
  return output(step, time, file_name, 2, (double []){div_max, div_sum});
}

static int monitor_quadratic_quantities(
    const size_t step,
    const double time,
    const domain_t * const domain,
    const flow_field_t * const flow_field
) {
  const char file_name[] = ROOT_DIRECTORY "quadratic_quantities.dat";
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const coordinate_t * const y_coordinate = &domain->y_coordinate;
  const double * const x_face_scale_factors = x_coordinate->face_scale_factors;
  const double * const x_center_scale_factors = x_coordinate->center_scale_factors;
  const double * const y_face_scale_factors = y_coordinate->face_scale_factors;
  const double * const y_center_scale_factors = y_coordinate->center_scale_factors;
  const array_t * const ux_array = flow_field->ux;
  const array_t * const uy_array = flow_field->uy;
  const array_t * const temperature_array = flow_field->temperature;
  const array_t * const salinity_array = flow_field->salinity;
  const double (* const ux)[NX + 2] = ux_array->buffer;
  const double (* const uy)[NX + 2] = uy_array->buffer;
  const double (* const temperature)[NX + 2] = temperature_array->buffer;
  const double (* const salinity)[NX + 2] = salinity_array->buffer;
  double quantities[3] = {0.};
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 2; i <= NX; i++) {
      const double hx = x_face_scale_factors[i];
      const double hy = y_center_scale_factors[j];
      quantities[0] += 0.5 * hx * hy * pow(ux[j][i], 2.);
    }
  }
  for (size_t j = 2; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double hx = x_center_scale_factors[i];
      const double hy = y_face_scale_factors[j];
      quantities[0] += 0.5 * hx * hy * pow(uy[j][i], 2.);
    }
  }
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double hx = x_center_scale_factors[i];
      const double hy = y_center_scale_factors[j];
      quantities[1] += 0.5 * hx * hy * pow(temperature[j][i], 2.);
    }
  }
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double hx = x_center_scale_factors[i];
      const double hy = y_center_scale_factors[j];
      quantities[2] += 0.5 * hx * hy * pow(salinity[j][i], 2.);
    }
  }
  return output(step, time, file_name, 3, quantities);
}

int monitor(
    const size_t step,
    const double time,
    const double dt,
    const domain_t * const domain,
    const flow_field_t * const flow_field
) {
  fprintf(
      stdout,
      "step %10zu time % .2e dt % .2e\n",
      step, time, dt
  );
  if (0 != monitor_divergence(step, time, domain, flow_field)) {
    goto abort;
  }
  if (0 != monitor_quadratic_quantities(step, time, domain, flow_field)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

