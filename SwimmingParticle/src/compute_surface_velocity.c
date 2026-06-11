#include "./compute_surface_velocity.h"

int compute_surface_velocity(
    const domain_t * const domain,
    const double * const surface_concentration,
    double * const surface_velocity
) {
  const size_t ny = domain->ny;
  // differentiate in the azimuthal direction
  //   in the frequency domain by multiplying by I k:
  // I k c_k = I k (r_k + I i_k)
  //         = - k i_k + I k r_k
  // NOTE: i_k = - i_{n - k}
  surface_velocity[     0] = 0.;
  surface_velocity[ny / 2] = 0.;
  for (size_t k = 1; k < ny / 2; k++) {
    const double real = + surface_concentration[     k];
    const double imag = - surface_concentration[ny - k];
    surface_velocity[     k] = + (- imag * k);
    surface_velocity[ny - k] = - (+ real * k);
  }
  return 0;
}

