#include <math.h>
#include <stddef.h>
#include "./monitor.h"

static int compute_max_velocity(
    const vector_t * const velocity,
    double * const max_velocity
) {
  *max_velocity = 0.;
  for (int j = 1; j <= NY; j++) {
    for (int i = 1; i <= NX; i++) {
      const double * const vij = (*velocity)[j][i];
      const double local_velocity = sqrt(
          + vij[0] * vij[0]
          + vij[1] * vij[1]
      );
      *max_velocity = fmax(*max_velocity, local_velocity);
    }
  }
  return 0;
}

static int compute_nusselt(
    const double l,
    const double kappa,
    const vector_t * const velocity,
    const scalar_t * const temperature,
    double * const nusselt
) {
  double q_ref = 0.;
  double q = 0.;
  for (size_t i = 1; i <= NX; i++) {
    for (size_t j = 0; j <= NY; j++) {
      const double vm = (*velocity)[j][i][1];
      const double vp = (*velocity)[j + 1][i][1];
      const double v = 0.5 * vm + 0.5 * vp;
      const double tm = (*temperature)[j][i];
      const double tp = (*temperature)[j + 1][i];
      const double t = 0.5 * tm + 0.5 * tp;
      q_ref += kappa / l;
      q += v * t - kappa * (tp - tm);
    }
  }
  *nusselt = q / q_ref;
  return 0;
}

int monitor(
    const double l,
    const double kappa,
    const vector_t * const velocity,
    const scalar_t * const temperature,
    double * const max_velocity,
    double * const nusselt
) {
  if (0 != compute_max_velocity(velocity, max_velocity)) {
    return 1;
  }
  if (0 != compute_nusselt(l, kappa, velocity, temperature, nusselt)) {
    return 1;
  }
  return 0;
}

