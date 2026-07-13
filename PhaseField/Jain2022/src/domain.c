#include "domain.h"

int domain_init(
    domain_t * const domain
) {
  const double lx = 1.;
  const double ly = 2.;
  const size_t nx = 64;
  const size_t ny = 128;
  const double dx = lx / nx;
  const double dy = ly / ny;
  domain->lx = lx;
  domain->ly = ly;
  domain->nx = nx;
  domain->ny = ny;
  domain->dx = dx;
  domain->dy = dy;
  return 0;
}

