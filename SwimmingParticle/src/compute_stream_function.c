#include <math.h>
#include <stddef.h>
#include "./compute_stream_function.h"
#include "./exchange_halo.h"

int compute_stream_function(
    const domain_t * const domain,
    rdft_plan_t * const rdft_plan,
    const double * const surface_velocity,
    double ** const stream_function
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const xf = domain->xf;
  // from surface velocity (and other boundary conditions),
  //   compute stream function in the azimuthal frequency space
  //   for each radial position
#pragma omp parallel for
  for (size_t i = 1; i <= nx + 1; i++) {
    double * const stream_function_r = stream_function[i];
    // k = 0 and ny / 2: the surface velocity is zero and thus zero
    stream_function_r[         1] = 0.;
    stream_function_r[ny / 2 + 1] = 0.;
    // other azimuthal wavenumbers
    const double x = xf[i];
    for (size_t k = 1; k < ny / 2; k++) {
      const double pref = 0.5 * (1. - pow(x, 2.)) / pow(x, k);
      double * const real = &stream_function_r[     k + 1];
      double * const imag = &stream_function_r[ny - k + 1];
      *real = pref * surface_velocity[     k];
      *imag = pref * surface_velocity[ny - k];
    }
  }
  // transform to physical space
#pragma omp parallel for
  for (size_t i = 1; i <= nx + 1; i++) {
    rdft_exec_b(rdft_plan, &stream_function[i][1]);
  }
  exchange_halo(domain, stream_function);
  return 0;
}

