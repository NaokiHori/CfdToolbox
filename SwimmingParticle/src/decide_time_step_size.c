#include <float.h>
#include <stddef.h>
#include <math.h>
#include "domain.h"
#include "time_marcher.h"
#include "./decide_time_step_size.h"

int decide_time_step_size(
    const double peclet,
    const domain_t * const domain,
    double ** const stream_function,
    double * const dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const hxxf = domain->hxxf;
  const double * const hxxc = domain->hxxc;
  const double * const hyxf = domain->hyxf;
  const double * const hyxc = domain->hyxc;
  const double small = 1e-8;
  double dt_adv = DBL_MAX;
  double dt_dif = DBL_MAX;
  *dt = 5e-1;
#pragma omp parallel for reduction(min: dt_adv)
  for (size_t i = 2; i <= nx; i++) {
    for (size_t j = 1; j <= ny; j++) {
      const double ux = + 1. / hyxf[i] * (
          - stream_function[i    ][j    ]
          + stream_function[i    ][j + 1]
      );
      dt_adv = fmin(dt_adv, hxxf[i] / fmax(small, fabs(ux)));
    }
  }
#pragma omp parallel for reduction(min: dt_adv)
  for (size_t i = 1; i <= nx; i++) {
    for (size_t j = 1; j <= ny; j++) {
      const double uy = - 1. / hxxc[i] * (
          - stream_function[i    ][j    ]
          + stream_function[i + 1][j    ]
      );
      dt_adv = fmin(dt_adv, hyxc[i] / fmax(small, fabs(uy)));
    }
  }
  if (0 == IMPLICIT_FACTOR) {
#pragma omp parallel for reduction(min: dt_dif)
    for (size_t i = 1; i <= nx; i++) {
      dt_dif = fmin(dt_dif, peclet * 0.5 / NDIMS * pow(hxxc[i], 2.));
      dt_dif = fmin(dt_dif, peclet * 0.5 / NDIMS * pow(hyxc[i], 2.));
    }
  } else {
    dt_dif = 1. / small;
  }
  const double safety_factor = 0.95;
  *dt = fmin(*dt, safety_factor * fmin(dt_adv, dt_dif));
  return 0;
}

