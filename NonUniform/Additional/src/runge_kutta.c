#include "runge_kutta.h"

const size_t N_RUNGE_KUTTA_STEP = 3;

int runge_kutta_coefficients_setup(
    const size_t step,
    runge_kutta_coefficients_t * const coefficients
) {
  if (0 == step) {
    coefficients->explicit_previous = 0. / 60.;
    coefficients->explicit_current = 32. / 60.;
  } else if (1 == step) {
    coefficients->explicit_previous = - 17. / 60.;
    coefficients->explicit_current = 25. / 60.;
  } else {
    coefficients->explicit_previous = - 25. / 60.;
    coefficients->explicit_current = 45. / 60.;
  }
  coefficients->implicit
    = coefficients->explicit_previous
    + coefficients->explicit_current;
  return 0;
}

