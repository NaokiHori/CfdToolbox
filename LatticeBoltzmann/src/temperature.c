#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "./halo.h"
#include "./output.h"
#include "./stream.h"
#include "./temperature.h"

static const double temperature_bottom = 0.5;
static const double temperature_top = -0.5;

static inline double compute_equilibrium(
    const size_t n,
    const double * const velocity,
    const double temperature
) {
  const double first =
    + 1. / CS2 * (double)c[n][0] * velocity[0]
    + 1. / CS2 * (double)c[n][1] * velocity[1];
  const double series = 1. + first;
  return w[n] * temperature * series;
}

static int impose_boundary_condition(
    const double tw,
    const double * const g0ij,
    const double * const g1ij,
    double * const g2ij,
    const double * const g3ij,
    const double * const g4ij,
    double * const g5ij,
    double * const g6ij,
    const double * const g7ij,
    const double * const g8ij
) {
  const double tp = 1. / (W1 + 2. * W2) * (
      tw
      - *g0ij - *g1ij - *g3ij
      - *g4ij - *g7ij - *g8ij
  );
  const double up = - 0.5 * CS2 / W2 * (
      + *g1ij - *g3ij
      - *g7ij + *g8ij
  ) / tp;
  *g2ij = compute_equilibrium(2, (double [NDIMS]){up, 0.}, tp);
  *g5ij = compute_equilibrium(5, (double [NDIMS]){up, 0.}, tp);
  *g6ij = compute_equilibrium(6, (double [NDIMS]){up, 0.}, tp);
  return 0;
}

static int impose_boundary_conditions(
    df_t * const g
) {
  // at y-negative wall
  for (size_t i = 1; i <= NX; i++) {
    double * const gij = (*g)[0][i];
    if (0 != impose_boundary_condition(
        temperature_bottom,
        gij + 0,
        gij + 1,
        gij + 2,
        gij + 3,
        gij + 4,
        gij + 5,
        gij + 6,
        gij + 7,
        gij + 8
    )) {
      return 1;
    }
  }
  // at y-positive wall
  for (size_t i = 1; i <= NX; i++) {
    double * const gij = (*g)[NY + 1][i];
    if (0 != impose_boundary_condition(
        temperature_top,
        gij + 0,
        gij + 3,
        gij + 4,
        gij + 1,
        gij + 2,
        gij + 7,
        gij + 8,
        gij + 5,
        gij + 6
    )) {
      return 1;
    }
  }
  return 0;
}

int initialize_temperature_distribution_function(
    df_t * const g
) {
  const double amp = 1e-4;
  const double vij[NDIMS] = {0., 0.};
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double r = 1. * rand() / RAND_MAX;
      const double tij = amp * ((temperature_top - temperature_bottom) * r + temperature_bottom);
      double * const gij = (*g)[j][i];
      for (size_t n = 0; n < NDIRS; n++) {
        gij[n] = compute_equilibrium(n, vij, tij);
      }
    }
  }
  for (size_t i = 1; i <= NX; i++) {
    double * const gij = (*g)[0][i];
    for (size_t n = 0; n < NDIRS; n++) {
      gij[n] = compute_equilibrium(n, vij, temperature_bottom);
    }
  }
  for (size_t i = 1; i <= NX; i++) {
    double * const gij = (*g)[NY + 1][i];
    for (size_t n = 0; n < NDIRS; n++) {
      gij[n] = compute_equilibrium(n, vij, temperature_top);
    }
  }
  if (0 != exchange_halo(g)) {
    return 1;
  }
  return 0;
}

int compute_temperature_macroscopic_field(
    const df_t * const g,
    scalar_t * const temperature
) {
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double * const gij = (*g)[j][i];
      double * const tij = &(*temperature)[j][i];
      *tij = 0.;
      for (size_t n = 0; n < NDIRS; n++) {
        *tij += gij[n];
      }
    }
  }
  return 0;
}

int process_temperature_collision(
    const double diffusivity,
    const vector_t * const velocity,
    const scalar_t * const temperature,
    df_t * const g
) {
  const double t_inv = 1. / (0.5 + diffusivity / CS2);
#pragma omp parallel for
  for (size_t j = 1; j <= NY; j++) {
    for (size_t i = 1; i <= NX; i++) {
      const double * const vij = (*velocity)[j][i];
      const double tij = (*temperature)[j][i];
      double * const gij = (*g)[j][i];
      for (size_t n = 0; n < NDIRS; n++) {
        double * const gnij = gij + n;
        const double g_eq = compute_equilibrium(n, vij, tij);
        *gnij = *gnij - t_inv * (*gnij - g_eq);
      }
    }
  }
  if (0 != exchange_halo(g)) {
    return 1;
  }
  return 0;
}

int process_temperature_streaming(
    df_t * const g
) {
  if (0 != process_streaming_inplace(g)) {
    return 1;
  }
  if (0 != impose_boundary_conditions(g)) {
    return 1;
  }
  if (0 != exchange_halo(g)) {
    return 1;
  }
  return 0;
}

