#include <math.h>
#include "array.h"
#include "boundary_condition.h"
#include "exchange_halo.h"
#include "interface_field.h"
#include "logger.h"

// diffusivity
const double interface_gamma = 1.;
// sharpness parameter
const double interface_epsilon = 1.;

static int init_vof(
    const domain_t * const domain,
    double ** const vof
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double lx = domain->lx;
  const double ly = domain->ly;
  const double epsilon = interface_epsilon * fmax(dx, dy);
  const double center_x = 0.5 * lx;
  const double center_y = 0.5 * ly;
  const double radius = 0.25 * fmin(lx, ly);
  for (size_t j = 1; j <= ny; j++) {
    const double y = 0.5 * (2 * j - 1) * dy;
    for (size_t i = 1; i <= nx; i++) {
      const double x = 0.5 * (2 * i - 1) * dx;
      const double distance = radius - sqrt(pow(x - center_x, 2.) + pow(y - center_y, 2.));
      vof[j][i] = 0.5 * (1. + tanh(0.5 / epsilon * distance));
    }
  }
  if (0 != impose_boundary_condition_vof_x(domain, vof)) {
    LOGGER_FAILURE("failed to impose boundary condition in x (vof)");
    goto abort;
  }
  if (0 != exchange_halo_y(domain, vof)) {
    LOGGER_FAILURE("failed to exchange halo in y (vof)");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

int interface_field_init(
    const domain_t * const domain,
    interface_field_t * const interface_field
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  array_init(nx + 2, ny + 2, &interface_field->vof);
  if (0 != init_vof(domain, interface_field->vof)) {
    LOGGER_FAILURE("failed to initialize vof");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

int interface_field_finalize(
    interface_field_t * const interface_field
) {
  array_finalize(&interface_field->vof);
  return 0;
}

