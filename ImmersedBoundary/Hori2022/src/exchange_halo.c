#include <stddef.h> // size_t
#include "logger.h" // LOGGER_FAILURE
#include "domain.h" // nx, ny
#include "exchange_halo.h"

int exchange_halo_y(
    const domain_t * const domain,
    double ** const array
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  for (size_t i = 0; i <= nx + 1; i++) {
    array[     0][i] = array[ny][i];
    array[ny + 1][i] = array[ 1][i];
  }
  return 0;
}

