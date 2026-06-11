#include <math.h>
#include <stddef.h>
#include "domain.h"
#include "./fluid.h"
#include "./halo.h"
#include "./output.h"
#include "./stream.h"

static inline double compute_equilibrium(
    const size_t n,
    const double density,
    const double * const velocity
) {
  const double first = 1. / CS2 * (
      + c[n][0] * velocity[0]
      + c[n][1] * velocity[1]
  );
  const double second =
    + 0.5 * first * first;
  const double third = - 0.5 / CS2 * (
      + velocity[0] * velocity[0]
      + velocity[1] * velocity[1]
  );
  const double series = 1. + first + second + third;
  return density * w[n] * series;
}

static int impose_boundary_condition(
    const double * const f1ij,
    double * const f2ij,
    const double * const f3ij,
    const double * const f4ij,
    double * const f5ij,
    double * const f6ij,
    const double * const f7ij,
    const double * const f8ij
) {
  const double rp = 1. / (W1 + 2. * W2) * (
      + *f4ij + *f7ij + *f8ij
  );
  const double up = 0.5 * CS2 / W2 / rp * (
      - *f1ij + *f3ij
      + *f7ij - *f8ij
  );
  *f2ij = compute_equilibrium(2, rp, (double [NDIMS]){up, 0.});
  *f5ij = compute_equilibrium(5, rp, (double [NDIMS]){up, 0.});
  *f6ij = compute_equilibrium(6, rp, (double [NDIMS]){up, 0.});
  return 0;
}

static int impose_boundary_conditions(
    df_t * const f
) {
  // at y-negative wall
  for (size_t i = 1; i <= NX; i++) {
    double * const fij = (*f)[0][i];
    if (0 != impose_boundary_condition(
        fij + 1,
        fij + 2,
        fij + 3,
        fij + 4,
        fij + 5,
        fij + 6,
        fij + 7,
        fij + 8
    )) {
      return 1;
    }
  }
  // at y-positive wall
  for (size_t i = 1; i <= NX; i++) {
    double * const fij = (*f)[NY + 1][i];
    if (0 != impose_boundary_condition(
        fij + 3,
        fij + 4,
        fij + 1,
        fij + 2,
        fij + 7,
        fij + 8,
        fij + 5,
        fij + 6
    )) {
      return 1;
    }
  }
  return 0;
}

int initialize_fluid_distribution_function(
    df_t * const f
) {
  const double dij = 1.;
  const double vij[NDIMS] = {0., 0.};
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      double * const fij = (*f)[j][i];
      for (size_t n = 0; n < NDIRS; n++) {
        fij[n] = compute_equilibrium(n, dij, vij);
      }
    }
  }
  for (size_t i = 1; i <= NX; i++) {
    double * const fij = (*f)[0][i];
    for (size_t n = 0; n < NDIRS; n++) {
      fij[n] = compute_equilibrium(n, dij, vij);
    }
  }
  for (size_t i = 1; i <= NX; i++) {
    double * const fij = (*f)[NY + 1][i];
    for (size_t n = 0; n < NDIRS; n++) {
      fij[n] = compute_equilibrium(n, dij, vij);
    }
  }
  if (0 != exchange_halo(f)) {
    return 1;
  }
  return 0;
}

int compute_fluid_macroscopic_field(
    const df_t * const f,
    scalar_t * const density,
    vector_t * const velocity
) {
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double * const fij = (*f)[j][i];
      double * const dij = &(*density)[j][i];
      double * const vij = (*velocity)[j][i];
      *dij = 0.;
      vij[0] = 0.;
      vij[1] = 0.;
      for (size_t n = 0; n < NDIRS; n++) {
        const double fnij = fij[n];
        const int * const cn = c[n];
        *dij += fnij;
        vij[0] += fnij * cn[0];
        vij[1] += fnij * cn[1];
      }
      // N.B. ignoring zero division
      vij[0] /= *dij;
      vij[1] /= *dij;
    }
  }
  return 0;
}

int process_fluid_collision(
    const double diffusivity,
    const double acceleration,
    const scalar_t * const density,
    const vector_t * const velocity,
    const scalar_t * const temperature,
    df_t * const f
) {
  const double t_inv = 1. / (0.5 + diffusivity / CS2);
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double dij = (*density)[j][i];
      const double * const vij = (*velocity)[j][i];
      const double tij = (*temperature)[j][i];
      const double ay = acceleration * tij;
      const double fy = dij * ay;
      // shifted velocity
      const double uij[NDIMS] = {
        vij[0],
        vij[1] + 0.5 * ay,
      };
      double * const fij = (*f)[j][i];
      for (size_t n = 0; n < NDIRS; n++) {
        double * const fnij = fij + n;
        const int * const cn = c[n];
        // collision
        const double f_eq = compute_equilibrium(n, dij, uij);
        *fnij = *fnij - t_inv * (*fnij - f_eq);
        // forcing
        *fnij += (1. - 0.5 * t_inv) / CS2 * w[n] * (
            + (cn[1] - uij[1])
            + 1. / CS2 * (
                + cn[0] * uij[0]
                + cn[1] * uij[1]
            ) * cn[1]
        ) * fy;
      }
    }
  }
  if (0 != exchange_halo(f)) {
    return 1;
  }
  return 0;
}

int process_fluid_streaming(
    df_t * const f
) {
  if (0 != process_streaming_inplace(f)) {
    return 1;
  }
  if (0 != impose_boundary_conditions(f)) {
    return 1;
  }
  if (0 != exchange_halo(f)) {
    return 1;
  }
  return 0;
}

