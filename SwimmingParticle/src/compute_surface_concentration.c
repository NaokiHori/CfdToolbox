#include <math.h>
#include <stddef.h>
#include "constant.h"
#include "./compute_surface_concentration.h"

int compute_surface_concentration(
    const domain_t * const domain,
    rdft_plan_t * const rdft_plan,
    double ** const concentration,
    double * const surface_concentration
) {
  // compute concentration on the object surface in the frequency domain
  //   using c_1, c_2, ..., c_{ny}
  const size_t ny = domain->ny;
  for (size_t j = 0; j < ny; j++) {
    surface_concentration[j] = concentration[0][j + 1];
  }
  rdft_exec_f(rdft_plan, surface_concentration);
  // since the concentration field is defined at each azimuthal cell center,
  //   we need to perform a phase shift by half grid size
  for (size_t k = 1; k < ny / 2; k++) {
    //   c_k exp(- pi k / n I)
    // = (r_k + I i_k) * (cos + I sin)
    // = (r_k * cos - i_k * sin) + I (i_k * cos + r_k * sin)
    // NOTE: i_k = - i_{n - k}
    const double real = + surface_concentration[     k];
    const double imag = - surface_concentration[ny - k];
    const double cosk = cos(- PI * k / ny);
    const double sink = sin(- PI * k / ny);
    surface_concentration[     k] = + (real * cosk - imag * sink);
    surface_concentration[ny - k] = - (imag * cosk + real * sink);
  }
  return 0;
}

