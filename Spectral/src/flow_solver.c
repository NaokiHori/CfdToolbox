#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "derivative.h"
#include "flow_solver.h"
#include "parameter.h"

int initialize_flow_solver(
    const domain_t * const domain,
    flow_solver_t * const flow_solver
) {
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  const double * const kxs = domain->kxs;
  const double * const kys = domain->kys;
  double ** const integrating_factor_coefficients = &flow_solver->integrating_factor_coefficients;
  *integrating_factor_coefficients = malloc(knx * kny * sizeof(double));
  for (size_t n = 0; n < knx * kny; n++) {
    const size_t i = n / kny;
    const size_t j = n % kny;
    const double kx = kxs[i];
    const double ky = kys[j];
    const double k_square = kx * kx + ky * ky;
    (*integrating_factor_coefficients)[n] = 1. / Re * k_square;
  }
  if (0 != initialize_runge_kutta(
      knx * kny,
      *integrating_factor_coefficients,
      &flow_solver->runge_kutta_ux
  )) {
    return 1;
  }
  if (0 != initialize_runge_kutta(
      knx * kny,
      *integrating_factor_coefficients,
      &flow_solver->runge_kutta_uy
  )) {
    return 1;
  }
  for (size_t n = 0; n < N_BUFFERS; n++) {
    flow_solver->buffers[n] = malloc(knx * kny * sizeof(double complex));
  }
  return 0;
}

int finalize_flow_solver(
    flow_solver_t * const flow_solver
) {
  free(flow_solver->integrating_factor_coefficients);
  if (0 != finalize_runge_kutta(&flow_solver->runge_kutta_ux)) {
    return 1;
  }
  if (0 != finalize_runge_kutta(&flow_solver->runge_kutta_uy)) {
    return 1;
  }
  for (size_t n = 0; n < N_BUFFERS; n++) {
    free(flow_solver->buffers[n]);
  }
  return 0;
}

static int dealias(
    const domain_t * const domain,
    const double complex * const original,
    double complex * const clipped
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
#pragma omp parallel for
  for (size_t n = 0; n < knx * kny; n++) {
    const size_t i = n / kny;
    const size_t j = n % kny;
    const size_t mode_j = (j < ny / 2 + 1) ? j : (ny - j);
    if (i < nx / 3 && mode_j < ny / 3) {
      clipped[n] = original[n];
    } else {
      clipped[n] = 0.;
    }
  }
  return 0;
}

static int compute_time_step_size(
    const domain_t * const domain,
    const double * const phys_ux,
    const double * const phys_uy,
    double * const dt
) {
  const double lx = domain->lx;
  const double ly = domain->ly;
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  double max_ux = 0.;
  double max_uy = 0.;
#pragma omp parallel for reduction(max: max_ux)
  for (size_t n = 0; n < nx * ny; n++) {
    const double ux = phys_ux[n];
    const double abs_ux = fabs(ux);
    if (max_ux < abs_ux) {
      max_ux = abs_ux;
    }
  }
#pragma omp parallel for reduction(max: max_uy)
  for (size_t n = 0; n < nx * ny; n++) {
    const double uy = phys_uy[n];
    const double abs_uy = fabs(uy);
    if (max_uy < abs_uy) {
      max_uy = abs_uy;
    }
  }
  const double safety_factor = 2.0e-1;
  const double grid_size = fmin(lx / nx, ly / ny);
  const double max_velocity = fmax(max_ux, max_uy);
  *dt = safety_factor * grid_size / max_velocity;
  return 0;
}

int integrate(
    const domain_t * const domain,
    transformer_t * const transformer,
    flow_solver_t * const flow_solver,
    flow_field_t * const flow_field,
    double * const dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  const double * const kxs = domain->kxs;
  const double * const kys = domain->kys;
  // main arrays: velocity in spectral domain
  double complex * const spec_ux = flow_field->ux;
  double complex * const spec_uy = flow_field->uy;
  // inverse transformed arrays
  double * const phys_ux = flow_solver->buffers[0];
  double * const phys_uy = flow_solver->buffers[1];
  double complex * const clipped_spec_ux = flow_solver->buffers[0];
  double complex * const clipped_spec_uy = flow_solver->buffers[1];
  // non-linear terms
  double * const phys_uxux = flow_solver->buffers[0];
  double * const phys_uxuy = flow_solver->buffers[2];
  double * const phys_uyuy = flow_solver->buffers[1];
  double complex * const spec_uxux = flow_solver->buffers[0];
  double complex * const spec_uxuy = flow_solver->buffers[2];
  double complex * const spec_uyuy = flow_solver->buffers[1];
  // advection terms (x components)
  double complex * const spec_duxuxdx = flow_solver->buffers[0];
  double complex * const spec_duyuxdy = flow_solver->buffers[3];
  // advection terms (y components)
  double complex * const spec_duxuydx = flow_solver->buffers[2];
  double complex * const spec_duyuydy = flow_solver->buffers[1];
  // runge-kutta integrators
  runge_kutta_t * const runge_kutta_ux = flow_solver->runge_kutta_ux;
  runge_kutta_t * const runge_kutta_uy = flow_solver->runge_kutta_uy;
  double complex * const right_hand_side_ux = flow_solver->buffers[0];
  double complex * const right_hand_side_uy = flow_solver->buffers[1];
  for (size_t runge_kutta_step = 0; runge_kutta_step < runge_kutta_step_max; runge_kutta_step++) {
    // go to physical space
    if (0 != dealias(domain, spec_ux, clipped_spec_ux)) {
      return 1;
    }
    if (0 != dealias(domain, spec_uy, clipped_spec_uy)) {
      return 1;
    }
    if (0 != transform_spectral_to_physical(transformer, clipped_spec_ux, phys_ux)) {
      return 1;
    }
    if (0 != transform_spectral_to_physical(transformer, clipped_spec_uy, phys_uy)) {
      return 1;
    }
    // set time-step size based on cfl condition
    if (0 == runge_kutta_step) {
      if (0 != compute_time_step_size(domain, phys_ux, phys_uy, dt)) {
        return 1;
      }
      if (0 != runge_kutta_set_time_step_size(runge_kutta_ux, *dt)) {
        return 1;
      }
      if (0 != runge_kutta_set_time_step_size(runge_kutta_uy, *dt)) {
        return 1;
      }
    }
    // evaluate non-linear terms
    // NOTE: order matters as buffers are reused
#pragma omp parallel for
    for (size_t n = 0; n < nx * ny; n++) {
      phys_uxuy[n] = phys_ux[n] * phys_uy[n];
    }
#pragma omp parallel for
    for (size_t n = 0; n < nx * ny; n++) {
      phys_uxux[n] = phys_ux[n] * phys_ux[n];
    }
#pragma omp parallel for
    for (size_t n = 0; n < nx * ny; n++) {
      phys_uyuy[n] = phys_uy[n] * phys_uy[n];
    }
    // go back to spectral domain
    if (0 != transform_physical_to_spectral(transformer, phys_uxux, spec_uxux)) {
      return 1;
    }
    if (0 != transform_physical_to_spectral(transformer, phys_uxuy, spec_uxuy)) {
      return 1;
    }
    if (0 != transform_physical_to_spectral(transformer, phys_uyuy, spec_uyuy)) {
      return 1;
    }
    // differentiate
    // NOTE: order matters as buffers are reused
    if (0 != compute_y_derivative(domain, 1., spec_uxuy, spec_duyuxdy)) {
      return 1;
    }
    if (0 != compute_x_derivative(domain, 1., spec_uxuy, spec_duxuydx)) {
      return 1;
    }
    if (0 != compute_x_derivative(domain, 1., spec_uxux, spec_duxuxdx)) {
      return 1;
    }
    if (0 != compute_y_derivative(domain, 1., spec_uyuy, spec_duyuydy)) {
      return 1;
    }
    // update spectral velocity
#pragma omp parallel for
    for (size_t n = 0; n < knx * kny; n++) {
      const size_t i = n / kny;
      const size_t j = n % kny;
      const double kx = kxs[i];
      const double ky = kys[j];
      // i = 0 && j = 0 will be modified afterwards
      const double k_square_inverse = 0 == i && 0 == j ? 0. : 1. / (kx * kx + ky * ky);
      const double complex advx = - spec_duxuxdx[n] - spec_duyuxdy[n];
      const double complex advy = - spec_duxuydx[n] - spec_duyuydy[n];
      const double complex k_adv = kx * advx + ky * advy;
      right_hand_side_ux[n] = advx - kx * k_square_inverse * k_adv;
      right_hand_side_uy[n] = advy - ky * k_square_inverse * k_adv;
    }
    if (0 != runge_kutta_update(runge_kutta_ux, right_hand_side_ux, spec_ux)) {
      return 1;
    }
    if (0 != runge_kutta_update(runge_kutta_uy, right_hand_side_uy, spec_uy)) {
      return 1;
    }
    // no net movement
    spec_ux[0] = 0.;
    spec_uy[0] = 0.;
    // nullify components at Nyquist frequencies
    for (size_t i = 0; i < knx; i++) {
      const size_t j = ny / 2;
      spec_ux[i * ny + j] = 0.;
      spec_uy[i * ny + j] = 0.;
    }
    for (size_t j = 0; j < kny; j++) {
      const size_t i = nx / 2;
      spec_ux[i * ny + j] = 0.;
      spec_uy[i * ny + j] = 0.;
    }
  }
  return 0;
}

