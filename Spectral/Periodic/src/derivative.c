#include <stddef.h>
#include "derivative.h"

int compute_x_derivative(
    const domain_t * const domain,
    const double sign,
    const double complex * const src,
    double complex * const dst
) {
  // multiply by i * kx
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  const double * const kxs = domain->kxs;
#pragma omp parallel for
  for (size_t i = 0; i < knx - 1; i++) {
    const double kx = kxs[i];
    for (size_t j = 0; j < kny; j++) {
      dst[i * kny + j] = sign * kx * I * src[i * kny + j];
    }
  }
  {
    const size_t i = knx - 1;
    for (size_t j = 0; j < kny; j++) {
      dst[i * kny + j] = 0.;
    }
  }
  return 0;
}

int compute_y_derivative(
    const domain_t * const domain,
    const double sign,
    const double complex * const src,
    double complex * const dst
) {
  // multiply by i * ky
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  const double * const kys = domain->kys;
#pragma omp parallel for
  for (size_t i = 0; i < knx; i++) {
    for (size_t j = 0; j < kny; j++) {
      if (kny / 2 == j) {
        dst[i * kny + j] = 0.;
      } else {
        const double ky = kys[j];
        dst[i * kny + j] = sign * ky * I * src[i * kny + j];
      }
    }
  }
  return 0;
}

#if defined(TEST_DERIVATIVE)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constant.h"

static int test_compute_x_derivative(
    void
) {
  domain_t domain = {0};
  initialize_domain(
      2. * PI, 2. * PI,
      4, 3,
      &domain
  );
  const size_t knx = domain.knx;
  const size_t kny = domain.kny;
  const double complex src[9] = {
    1100. - 11. * I, 2100. - 21. * I, 3100. - 31. * I,
    1200. - 12. * I, 2200. - 22. * I, 3200. - 32. * I,
    1300. - 13. * I, 2300. - 23. * I, 3300. - 33. * I,
  };
  double complex dst[9] = {0.};
  // wavenumbers are multiplied
  // Nyquist frequencies are nullified
  const double complex expected[] = {
     0. +    0. * I,  0. +    0. * I,  0. +    0. * I,
    12. + 1200. * I, 22. + 2200. * I, 32. + 3200. * I,
     0. +    0. * I,  0. +    0. * I,  0. +    0. * I,
  };
  compute_x_derivative(&domain, 1., src, dst);
  double error = 0.;
  for (size_t n = 0; n < knx * kny; n++) {
    error += pow(cabs(dst[n] - expected[n]), 2.);
  }
  error = sqrt(error / knx / kny);
  printf("% .1e\n", error);
  finalize_domain(&domain);
  return 0;
}

static int test_compute_y_derivative(
    void
) {
  domain_t domain = {0};
  initialize_domain(
      2. * PI, 2. * PI,
      4, 4,
      &domain
  );
  const size_t knx = domain.knx;
  const size_t kny = domain.kny;
  const double complex src[12] = {
    1100. - 11. * I, 2100. - 21. * I, 3100. - 31. * I, 4100. - 41. * I,
    1200. - 12. * I, 2200. - 22. * I, 3200. - 32. * I, 4200. - 42. * I,
    1300. - 13. * I, 2300. - 23. * I, 3300. - 33. * I, 4300. - 43. * I,
  };
  double complex dst[12] = {0.};
  // wavenumbers are multiplied
  // Nyquist frequencies are nullified
  const double complex expected[12] = {
    0. + 0. * I, 21. + 2100. * I, 0. + 0. * I, - 41. - 4100. * I,
    0. + 0. * I, 22. + 2200. * I, 0. + 0. * I, - 42. - 4200. * I,
    0. + 0. * I, 23. + 2300. * I, 0. + 0. * I, - 43. - 4300. * I,
  };
  compute_y_derivative(&domain, 1., src, dst);
  double error = 0.;
  for (size_t n = 0; n < knx * kny; n++) {
    error += pow(cabs(dst[n] - expected[n]), 2.);
  }
  error = sqrt(error / knx / kny);
  printf("% .1e\n", error);
  finalize_domain(&domain);
  return 0;
}

int main(
    void
) {
  test_compute_x_derivative();
  test_compute_y_derivative();
  return 0;
}

#endif // TEST_DERIVATIVE

