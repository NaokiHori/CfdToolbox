#include <stdlib.h>
#include "array.h"
#include "boundary_condition.h"
#include "domain.h"
#include "flow_field.h"
#include "logger.h"
#include "memory.h"

static int init_ux(
    array_t * const ux_array
) {
  double (* const ux)[NX + 2] = ux_array->buffer;
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      ux[j][i] = 0.;
    }
  }
  if (0 != impose_boundary_condition_ux_x(ux_array)) {
    goto abort;
  }
  if (0 != impose_boundary_condition_ux_y(ux_array)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

static int init_uy(
    array_t * const uy_array
) {
  double (* const uy)[NX + 2] = uy_array->buffer;
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      uy[j][i] = 0.;
    }
  }
  if (0 != impose_boundary_condition_uy_x(uy_array)) {
    goto abort;
  }
  if (0 != impose_boundary_condition_uy_y(uy_array)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

static int init_p(
    array_t * const p_array
) {
  double (* const p)[NX + 2] = p_array->buffer;
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      p[j][i] = 0.;
    }
  }
  // boundary conditions are implied in Poisson equation
  return 0;
}

static int init_scalar(
    const domain_t * const domain,
    array_t * const scalar_array
) {
  const coordinate_t * const x_coordinate = &domain->x_coordinate;
  const double lx = x_coordinate->length;
  const double * const x_centers = x_coordinate->centers;
  double (* const scalar)[NX + 2] = scalar_array->buffer;
  // perturbed randomly
  const double magnitude = 1e-2;
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double x_center = x_centers[i];
      scalar[j][i] = magnitude * (- 0.5 + 1. * rand() / RAND_MAX - 0.5 + x_center / lx);
    }
  }
  if (0 != impose_boundary_condition_scalar_x(scalar_array)) {
    goto abort;
  }
  if (0 != impose_boundary_condition_scalar_y(scalar_array)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

int flow_field_init(
    const domain_t * const domain,
    flow_field_t * const flow_field
) {
  flow_field->ux = memory_alloc(1, sizeof(array_t));
  flow_field->uy = memory_alloc(1, sizeof(array_t));
  flow_field->p = memory_alloc(1, sizeof(array_t));
  flow_field->temperature = memory_alloc(1, sizeof(array_t));
  flow_field->salinity = memory_alloc(1, sizeof(array_t));
  if (0 != init_ux(flow_field->ux)) {
    goto abort;
  }
  if (0 != init_uy(flow_field->uy)) {
    goto abort;
  }
  if (0 != init_p(flow_field->p)) {
    goto abort;
  }
  if (0 != init_scalar(domain, flow_field->temperature)) {
    goto abort;
  }
  if (0 != init_scalar(domain, flow_field->salinity)) {
    goto abort;
  }
  return 0;
abort:
  LOGGER_FAILURE();
  return 1;
}

int flow_field_finalize(
    flow_field_t * const flow_field
) {
  memory_free(flow_field->ux);
  memory_free(flow_field->uy);
  memory_free(flow_field->p);
  memory_free(flow_field->temperature);
  memory_free(flow_field->salinity);
  return 0;
}

