#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include "array.h"
#include "constant.h"
#include "./exchange_halo.h"
#include "./impose_boundary_condition.h"
#include "./time_marcher.h"
#include "./update_concentration_field.h"

int update_concentration_field(
    const double peclet,
    const domain_t * const domain,
    rdft_plan_t * const rdft_plan,
    tri_diagonal_matrix_t * const tri_diagonal_matrix,
    const double dt,
    const size_t runge_kutta_step,
    runge_kutta_buffer_t * const runge_kutta_buffer,
    double ** concentration
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const hxxf = domain->hxxf;
  const double * const hyxc = domain->hyxc;
  const double * const jdxf = domain->jdxf;
  const double * const jdxc = domain->jdxc;
  double ** const non_linear = runge_kutta_buffer->non_linear;
  double ** const increment = runge_kutta_buffer->increment;
  const double rk_a = RUNGE_KUTTA_COEFS.alpha[runge_kutta_step];
  const double rk_g = RUNGE_KUTTA_COEFS.gamma[runge_kutta_step];
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t j = 1; j <= ny; j++) {
      increment[i][j] =
        + rk_g * dt * non_linear[i][j]
        + rk_a * dt * increment[i][j];
    }
  }
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    rdft_exec_f(rdft_plan, &increment[i][1]);
  }
  const double pref = IMPLICIT_FACTOR * rk_a * dt / peclet;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    double * const l = tri_diagonal_matrix->l[j];
    double * const c = tri_diagonal_matrix->c[j];
    double * const u = tri_diagonal_matrix->u[j];
    for (size_t i = 1; i <= nx; i++) {
      l[i] = - pref / jdxc[i    ] * jdxf[i    ] / hxxf[i    ] / hxxf[i    ];
      u[i] = - pref / jdxc[i    ] * jdxf[i + 1] / hxxf[i + 1] / hxxf[i + 1];
      // radial Laplacian contribution
      // NOTE: impose Neumann boundary condition at particle surface
      c[i] = 1 == i ? - u[i] : - l[i] - u[i];
      // azimuthal Laplacian contribution
      c[i] += pref * pow(2. / hyxc[i    ] * sin(PI * (j - 1) / ny), 2.);
      // unity from the Helmholtz equation
      c[i] += 1.;
    }
    // forward sweep
    // NOTE: assuming non singular matrix (no zero division occurs)
    u[1] /= c[1];
    increment[1][j] /= c[1];
    for (size_t i = 2; i <= nx; i++) {
      const double val = 1. / (c[i] - l[i] * u[i - 1]);
      u[i] = val * u[i];
      increment[i][j] = val * (increment[i][j] - l[i] * increment[i - 1][j]);
    }
    // backward substitution
    for (size_t i = nx - 1; ; i--) {
      increment[i][j] -= u[i] * increment[i + 1][j];
      if (1 == i) {
        break;
      }
    }
  }
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    rdft_exec_b(rdft_plan, &increment[i][1]);
  }
#pragma omp parallel for
  for (size_t i = 1; i <= nx; i++) {
    for (size_t j = 1; j <= ny; j++) {
      concentration[i][j] += increment[i][j];
    }
  }
  impose_boundary_condition(domain, concentration);
  exchange_halo(domain, concentration);
  return 0;
}

