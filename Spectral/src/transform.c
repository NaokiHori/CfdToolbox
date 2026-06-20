#include <complex.h>
#include <stddef.h>
#include <stdlib.h>
#include "transform.h"
#include "./transform/dft.h"
#include "./transform/rdft.h"

#define N_BUFFERS 2

struct transformer_t {
  size_t nx;
  size_t ny;
  void * buffers[N_BUFFERS];
  rdft_plan_t * xplan;
  dft_plan_t * yplan;
};

int initialize_transformer(
    const size_t nx,
    const size_t ny,
    transformer_t ** const transformer
) {
  *transformer = malloc(sizeof(transformer_t));
  (*transformer)->nx = nx;
  (*transformer)->ny = ny;
  if (0 != rdft_init_plan(nx, &(*transformer)->xplan)) {
    return 1;
  }
  if (0 != dft_init_plan(ny, &(*transformer)->yplan)) {
    return 1;
  }
  const size_t buffer_size = (nx / 2 + 1) * ny * sizeof(double complex);
  for (size_t n = 0; n < N_BUFFERS; n++) {
    (*transformer)->buffers[n] = malloc(buffer_size);
  }
  return 0;
}

int finalize_transformer(
    transformer_t ** const transformer
) {
  rdft_destroy_plan(&(*transformer)->xplan);
  dft_destroy_plan(&(*transformer)->yplan);
  for (size_t n = 0; n < N_BUFFERS; n++) {
    free((*transformer)->buffers[n]);
  }
  free(*transformer);
  *transformer = NULL;
  return 0;
}

int transform_physical_to_spectral(
    transformer_t * const transformer,
    double * const physical,
    double complex * const spectral
) {
  const size_t nx = transformer->nx;
  const size_t ny = transformer->ny;
  rdft_plan_t * const xplan = transformer->xplan;
  dft_plan_t * const yplan = transformer->yplan;
  double * const buffer0 = transformer->buffers[0];
  double complex * const buffer1 = transformer->buffers[1];
  // x FFT
#pragma omp parallel for
  for (size_t j = 0; j < ny; j++) {
    rdft_exec_f(xplan, physical + j * nx, buffer0 + j * nx);
  }
  // transpose and pack as complex array
  {
    {
      const size_t i = 0;
      for (size_t j = 0; j < ny; j++) {
        buffer1[i * ny + j] = buffer0[j * nx + i];
      }
    }
#pragma omp parallel for
    for (size_t i = 1; i < nx / 2; i++) {
      for (size_t j = 0; j < ny; j++) {
        const double real = buffer0[j * nx + i];
        const double imag = - buffer0[j * nx + nx - i];
        buffer1[i * ny + j] = real + imag * I;
      }
    }
    {
      const size_t i = nx / 2;
      for (size_t j = 0; j < ny; j++) {
        buffer1[i * ny + j] = buffer0[j * nx + i];
      }
    }
  }
  // y FFT
#pragma omp parallel for
  for (size_t i = 0; i < nx / 2 + 1; i++) {
    dft_exec_f(yplan, buffer1 + i * ny, spectral + i * ny);
  }
  return 0;
}

int transform_spectral_to_physical(
    transformer_t * const transformer,
    const double complex * const spectral,
    double * const physical
) {
  const size_t nx = transformer->nx;
  const size_t ny = transformer->ny;
  rdft_plan_t * const xplan = transformer->xplan;
  dft_plan_t * const yplan = transformer->yplan;
  double complex * const buffer0 = transformer->buffers[0];
  double * const buffer1 = transformer->buffers[1];
  // y iFFT
#pragma omp parallel for
  for (size_t i = 0; i < nx / 2 + 1; i++) {
    dft_exec_b(yplan, spectral + i * ny, buffer0 + i * ny);
  }
  // transpose and pack as real array
  {
    {
      const size_t i = 0;
      for (size_t j = 0; j < ny; j++) {
        buffer1[j * nx + i] = creal(buffer0[i * ny + j]);
      }
    }
#pragma omp parallel for
    for (size_t i = 1; i < nx / 2; i++) {
      for (size_t j = 0; j < ny; j++) {
        const double complex value = buffer0[i * ny + j];
        buffer1[j * nx + i] = creal(value);
        buffer1[j * nx + nx - i] = - cimag(value);
      }
    }
    {
      const size_t i = nx / 2;
      for (size_t j = 0; j < ny; j++) {
        buffer1[j * nx + i] = creal(buffer0[i * ny + j]);
      }
    }
  }
  // x iFFT
#pragma omp parallel for
  for (size_t j = 0; j < ny; j++) {
    rdft_exec_b(xplan, buffer1 + j * nx, physical + j * nx);
  }
  return 0;
}

#if defined(TEST_TRANSFORM)

#include <math.h>
#include <stdio.h>
#include "constant.h"

static int test_transform_physical_to_spectral(
    void
) {
  const size_t nx = 4;
  const size_t ny = 4;
  transformer_t * transformer = NULL;
  initialize_transformer(nx, ny, &transformer);
  double * const physical = malloc(nx * ny * sizeof(double));
  double complex * const spectral = malloc((nx / 2 + 1) * ny * sizeof(double complex));
  for (size_t j = 0; j < ny; j++) {
    const double y = 2. * PI * j / ny;
    for (size_t i = 0; i < nx; i++) {
      const double x = 2. * PI * i / nx;
      physical[j * nx + i] = sin(x) * sin(y);
    }
  }
  const double complex expected[] = {
    + 0. + 0. * I, + 0. + 0. * I, + 0. + 0. * I, + 0. + 0. * I,
    + 0. + 0. * I, - 4. + 0. * I, + 0. + 0. * I, + 4. + 0. * I,
    + 0. + 0. * I, + 0. + 0. * I, + 0. + 0. * I, + 0. + 0. * I,
  };
  transform_physical_to_spectral(transformer, physical, spectral);
  double error = 0.;
  for (size_t n = 0; n < (nx / 2 + 1) * ny; n++) {
    error += pow(cabs(spectral[n] - expected[n]), 2.);
  }
  printf("% .1e\n", sqrt(error / (nx / 2 + 1) / ny));
  finalize_transformer(&transformer);
  free(physical);
  free(spectral);
  return 0;
}

static int test_transform_spectral_to_physical(
    void
) {
  const size_t nx = 4;
  const size_t ny = 4;
  transformer_t * transformer = NULL;
  initialize_transformer(nx, ny, &transformer);
  const double complex spectral[] = {
    + 0. + 0. * I, + 0. + 0. * I, + 0. + 0. * I, + 0. + 0. * I,
    + 0. + 0. * I, - 4. + 0. * I, + 0. + 0. * I, + 4. + 0. * I,
    + 0. + 0. * I, + 0. + 0. * I, + 0. + 0. * I, + 0. + 0. * I,
  };
  double * const physical = malloc(nx * ny * sizeof(double));
  double * const expected = malloc(nx * ny * sizeof(double));
  for (size_t j = 0; j < ny; j++) {
    const double y = 2. * PI * j / ny;
    for (size_t i = 0; i < nx; i++) {
      const double x = 2. * PI * i / nx;
      expected[j * nx + i] = sin(x) * sin(y);
    }
  }
  transform_spectral_to_physical(transformer, spectral, physical);
  double error = 0.;
  for (size_t n = 0; n < nx * ny; n++) {
    error += pow(physical[n] - expected[n], 2.);
  }
  printf("% .1e\n", sqrt(error / nx / ny));
  finalize_transformer(&transformer);
  free(physical);
  free(expected);
  return 0;
}

int main(
    void
) {
  test_transform_physical_to_spectral();
  test_transform_spectral_to_physical();
  return 0;
}

#endif // TEST_TRANSFORM

