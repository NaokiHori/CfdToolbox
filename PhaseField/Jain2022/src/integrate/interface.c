#include <math.h>
#include "boundary_condition.h"
#include "exchange_halo.h"
#include "logger.h"
#include "./interface.h"

int update_interface_field(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    interface_field_t * const interface_field,
    interface_solver_t * const interface_solver,
    const double dt
) {
  // a general-purpose small number to avoid singularities
  const double small = 1e-16;
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double dx = domain->dx;
  const double dy = domain->dy;
  double ** const ux = flow_field->ux;
  double ** const uy = flow_field->uy;
  double ** const vof = interface_field->vof;
  double ** const sdf = interface_solver->sdf;
  double ** const flux_x = interface_solver->flux_x;
  double ** const flux_y = interface_solver->flux_y;
  double ** const normal_x = interface_solver->normal_x;
  double ** const normal_y = interface_solver->normal_y;
  double ** const curvature = interface_solver->curvature;
  const double epsilon = interface_epsilon * fmax(dx, dy);
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 0; i <= nx + 1; i++) {
      const double local_vof = fmax(small, fmin(1. - small, vof[j][i]));
      sdf[j][i] = epsilon * log(local_vof / (1. - local_vof));
    }
  }
  if (0 != exchange_halo_y(domain, sdf)) {
    LOGGER_FAILURE("failed to exchange halo in y (sdf)");
    goto abort;
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double grad_x = (- sdf[j][i - 1] + sdf[j][i + 1]) / (2. * dx);
      const double grad_y = (- sdf[j - 1][i] + sdf[j + 1][i]) / (2. * dy);
      const double norm = fmax(small, sqrt(pow(grad_x, 2.) + pow(grad_y, 2.)));
      normal_x[j][i] = grad_x / norm;
      normal_y[j][i] = grad_y / norm;
    }
  }
  if (0 != impose_boundary_condition_normal_x_x(domain, normal_x)) {
    LOGGER_FAILURE("failed to impose boundary condition in x (normal_x)");
    goto abort;
  }
  if (0 != exchange_halo_y(domain, normal_x)) {
    LOGGER_FAILURE("failed to exchange halo in y (normal_x)");
    goto abort;
  }
  if (0 != exchange_halo_y(domain, normal_y)) {
    LOGGER_FAILURE("failed to exchange halo in y (normal_y)");
    goto abort;
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double div_normal_x = (- normal_x[j][i - 1] + normal_x[j][i + 1]) / (2. * dx);
      const double div_normal_y = (- normal_y[j - 1][i] + normal_y[j + 1][i]) / (2. * dy);
      curvature[j][i] = - div_normal_x - div_normal_y;
    }
  }
  if (0 != exchange_halo_y(domain, curvature)) {
    LOGGER_FAILURE("failed to exchange halo in y (curvature)");
    goto abort;
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    flux_x[j][1] = 0.;
    for (size_t i = 2; i <= nx; i++) {
      const double advection = ux[j][i] * (0.5 * vof[j][i - 1] + 0.5 * vof[j][i]);
      const double diffusion =
        - interface_gamma
          * epsilon
          * (- vof[j][i - 1] + vof[j][i]) / dx
        + 0.25
          * interface_gamma
          * (1. - pow(tanh(0.5 / epsilon * (0.5 * sdf[j][i - 1] + 0.5 * sdf[j][i])), 2.))
          * (0.5 * normal_x[j][i - 1] + 0.5 * normal_x[j][i]);
      flux_x[j][i] = advection + diffusion;
    }
    flux_x[j][nx + 1] = 0.;
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      const double advection = uy[j][i] * (0.5 * vof[j - 1][i] + 0.5 * vof[j][i]);
      const double diffusion =
        - interface_gamma
          * epsilon
          * (- vof[j - 1][i] + vof[j][i]) / dy
        + 0.25
          * interface_gamma
          * (1. - pow(tanh(0.5 / epsilon * (0.5 * sdf[j - 1][i] + 0.5 * sdf[j][i])), 2.))
          * (0.5 * normal_y[j - 1][i] + 0.5 * normal_y[j][i]);
      flux_y[j][i] = advection + diffusion;
    }
  }
  if (0 != exchange_halo_y(domain, flux_y)) {
    LOGGER_FAILURE("failed to exchange halo in y (flux_y)");
    goto abort;
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      vof[j][i] -= (- flux_x[j][i] + flux_x[j][i + 1]) / dx * dt;
      vof[j][i] -= (- flux_y[j][i] + flux_y[j + 1][i]) / dy * dt;
    }
  }
  if (0 != impose_boundary_condition_vof_x(domain, vof)) {
    LOGGER_FAILURE("failed to impose boundary condition in x (vof)");
    goto abort;
  }
  if (0 != exchange_halo_y(domain, vof)) {
    LOGGER_FAILURE("failed to exchange halo in y (vof)");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

