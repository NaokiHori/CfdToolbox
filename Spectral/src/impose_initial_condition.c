#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constant.h"
#include "derivative.h"
#include "output.h"
#include "./impose_initial_condition.h"

static double compute_gaussian_distribution(
    const double mean_x,
    const double mean_y,
    const double variance,
    const double x,
    const double y
) {
  const double prefactor = 1. / sqrt(2. * PI * variance);
  const double square_distance =
    + pow(x - mean_x, 2.)
    + pow(y - mean_y, 2.);
  return prefactor * exp(- 0.5 / variance * square_distance);
}

static int convert_streamfunction_to_velocity(
    const domain_t * const domain,
    double complex * const streamfunction,
    double complex * const ux,
    double complex * const uy
) {
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  const double * const kxs = domain->kxs;
  const double * const kys = domain->kys;
  for (size_t i = 0; i < knx; i++) {
    const double kx = kxs[i];
    for (size_t j = 0; j < kny; j++) {
      const double ky = kys[j];
      const double k_square_inverse = (0 == i && 0 == j) ? 0. : 1. / (kx * kx + ky * ky);
      streamfunction[i * kny + j] *= k_square_inverse;
    }
  }
  compute_y_derivative(domain, - 1., streamfunction, ux);
  compute_x_derivative(domain, + 1., streamfunction, uy);
  return 0;
}

int impose_initial_condition(
    const domain_t * const domain,
    transformer_t * const transformer,
    flow_field_t * const flow_field
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  double * const phys_streamfunction = malloc(nx * ny * sizeof(double));
  double complex * const spec_streamfunction = malloc(knx * kny * sizeof(double complex));
  double * const xs = malloc(nx * sizeof(double));
  double * const ys = malloc(ny * sizeof(double));
  for (size_t i = 0; i < nx; i++) {
    xs[i] = lx * i / nx - 0.5 * lx;
  }
  for (size_t j = 0; j < ny; j++) {
    ys[j] = ly * j / ny - 0.5 * ly;
  }
  output(OUTPUT_DIRECTORY "x.npy", 1, &nx, xs);
  output(OUTPUT_DIRECTORY "y.npy", 1, &ny, ys);
  const double l = fmin(lx, ly);
  for (size_t j = 0; j < ny; j++) {
    const double y = ys[j];
    for (size_t i = 0; i < nx; i++) {
      const double x = xs[i];
      // ref: https://basilisk.fr/src/test/stream.c
      phys_streamfunction[j * nx + i] =
        + 0.5 * compute_gaussian_distribution(- 0.1 * l, 0., 1e-3, x, y)
        + 0.5 * compute_gaussian_distribution(+ 0.1 * l, 0., 1e-3, x, y);
    }
  }
  free(xs);
  free(ys);
  if (0 != transform_physical_to_spectral(transformer, phys_streamfunction, spec_streamfunction)) {
    return 1;
  }
  if (0 != convert_streamfunction_to_velocity(
      domain,
      spec_streamfunction,
      flow_field->ux,
      flow_field->uy
  )) {
    return 1;
  }
  free(phys_streamfunction);
  free(spec_streamfunction);
  return 0;
}

