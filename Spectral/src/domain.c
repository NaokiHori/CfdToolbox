#include <stdlib.h>
#include "constant.h"
#include "domain.h"

int initialize_domain(
    const double lx,
    const double ly,
    const size_t nx,
    const size_t ny,
    domain_t * const domain
) {
  domain->lx = lx;
  domain->ly = ly;
  domain->nx = nx;
  domain->ny = ny;
  const size_t knx = nx / 2 + 1;
  const size_t kny = ny;
  domain->knx = knx;
  domain->kny = kny;
  double ** const kxs = &domain->kxs;
  double ** const kys = &domain->kys;
  *kxs = malloc(knx * sizeof(double));
  *kys = malloc(kny * sizeof(double));
  for (size_t i = 0; i < knx; i++) {
    (*kxs)[i] = 2. * PI / lx * i;
  }
  for (size_t j = 0; j < kny; j++) {
    (*kys)[j] = 2. * PI / ly * (j < ny / 2 + 1 ? (double)j : (double)j - (double)ny);
  }
  return 0;
}

int finalize_domain(
    domain_t * const domain
) {
  free(domain->kxs);
  free(domain->kys);
  return 0;
}

#if defined(TEST_DOMAIN)

#include <math.h>
#include <stdio.h>

static int test_kxs(
    const domain_t * const domain
) {
  const double expected[3] = {
    + 0. * 2. * PI,
    + 1. * 2. * PI,
    + 2. * 2. * PI,
  };
  double error = 0.;
  for (size_t i = 0; i < 3; i++) {
    error += pow(domain->kxs[i] - expected[i], 2.);
  }
  printf("% .1e\n", sqrt(error / 3));
  return 0;
}

static int test_kys(
    const domain_t * const domain
) {
  const double expected[4] = {
    + 0. * 2. * PI,
    + 1. * 2. * PI,
    + 2. * 2. * PI,
    - 1. * 2. * PI,
  };
  double error = 0.;
  for (size_t j = 0; j < 3; j++) {
    error += pow(domain->kys[j] - expected[j], 2.);
  }
  printf("% .1e\n", sqrt(error / 4));
  return 0;
}

int main(
    void
) {
  domain_t domain = {0};
  initialize_domain(1., 1., 4, 4, &domain);
  test_kxs(&domain);
  test_kys(&domain);
  finalize_domain(&domain);
  return 0;
}

#endif // TEST_DOMAIN

