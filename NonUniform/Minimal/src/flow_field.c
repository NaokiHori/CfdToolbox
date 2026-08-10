#include "logger.h"
#include "array.h"
#include "domain.h"
#include "flow_field.h"
#include "boundary_condition.h"

static int init_ux(
    const domain_t * const domain,
    double ** const ux
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      ux[j][i] = 0.;
    }
  }
  if (0 != impose_boundary_condition_ux_x(domain, ux)) {
    LOGGER_FAILURE("failed to impose boundary condition in x (ux)");
    goto abort;
  }
  if (0 != impose_boundary_condition_ux_y(domain, ux)) {
    LOGGER_FAILURE("failed to impose boundary condition in y (ux)");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

static int init_uy(
    const domain_t * const domain,
    double ** const uy
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      uy[j][i] = 0.;
    }
  }
  if (0 != impose_boundary_condition_uy_x(domain, uy)) {
    LOGGER_FAILURE("failed to impose boundary condition in x (uy)");
    goto abort;
  }
  if (0 != impose_boundary_condition_uy_y(domain, uy)) {
    LOGGER_FAILURE("failed to impose boundary condition in y (uy)");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

static int init_p(
    const domain_t * const domain,
    double ** const p
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      p[j][i] = 0.;
    }
  }
  return 0;
}

int flow_field_init(
    const domain_t * const domain,
    flow_field_t * const flow_field
) {
  const size_t nx = domain->x_coordinate.nitems;
  const size_t ny = domain->y_coordinate.nitems;
  array_init(nx + 2, ny + 2, &flow_field->ux);
  array_init(nx + 2, ny + 2, &flow_field->uy);
  array_init(nx + 2, ny + 2, &flow_field->p);
  if (0 != init_ux(domain, flow_field->ux)) {
    LOGGER_FAILURE("failed to initialize ux");
    goto abort;
  }
  if (0 != init_uy(domain, flow_field->uy)) {
    LOGGER_FAILURE("failed to initialize uy");
    goto abort;
  }
  if (0 != init_p(domain, flow_field->p)) {
    LOGGER_FAILURE("failed to initialize p");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

int flow_field_finalize(
    flow_field_t * const flow_field
) {
  array_finalize(&flow_field->ux);
  array_finalize(&flow_field->uy);
  array_finalize(&flow_field->p);
  return 0;
}

