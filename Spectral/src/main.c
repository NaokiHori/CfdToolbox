#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "constant.h"
#include "derivative.h"
#include "domain.h"
#include "flow_field.h"
#include "flow_solver.h"
#include "output.h"
#include "transform.h"
#include "./impose_initial_condition.h"

static int check_divergence(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    double * const max_divergence
) {
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  const double complex * const ux = flow_field->ux;
  const double complex * const uy = flow_field->uy;
  double complex * const duxdx = malloc(knx * kny * sizeof(double complex));
  double complex * const duydy = malloc(knx * kny * sizeof(double complex));
  compute_x_derivative(domain, 1., ux, duxdx);
  compute_y_derivative(domain, 1., uy, duydy);
  *max_divergence = 0;
  for (size_t n = 0; n < knx * kny; n++) {
    const double local_divergence = duxdx[n] + duydy[n];
    *max_divergence = fmax(*max_divergence, fabs(local_divergence));
  }
  free(duxdx);
  free(duydy);
  return 0;
}

static int output_vorticity(
    const char file_name[],
    const domain_t * const domain,
    transformer_t * const transformer,
    const flow_field_t * const flow_field
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  const double complex * const ux = flow_field->ux;
  const double complex * const uy = flow_field->uy;
  double * const phys_vorticity = malloc(nx * ny * sizeof(double));
  double complex * const spec_vorticity = malloc(knx * kny * sizeof(double complex));
  double complex * const spec_duydx = malloc(knx * kny * sizeof(double complex));
  double complex * const spec_duxdy = malloc(knx * kny * sizeof(double complex));
  compute_x_derivative(domain, 1., uy, spec_duydx);
  compute_y_derivative(domain, 1., ux, spec_duxdy);
  for (size_t n = 0; n < knx * kny; n++) {
    spec_vorticity[n] = spec_duydx[n] - spec_duxdy[n];
  }
  transform_spectral_to_physical(transformer, spec_vorticity, phys_vorticity);
  output(file_name, 2, (size_t []){ny, nx}, phys_vorticity);
  free(spec_duydx);
  free(spec_duxdy);
  free(spec_vorticity);
  free(phys_vorticity);
  return 0;
}

int main(
    void
) {
  const double lx = 1.;
  const double ly = 1.;
  const size_t nx = 256;
  const size_t ny = 256;
  domain_t domain = {0};
  transformer_t * transformer = NULL;
  flow_field_t flow_field = {0};
  flow_solver_t flow_solver = {0};
  if (0 != initialize_domain(lx, ly, nx, ny, &domain)) {
    return 1;
  }
  if (0 != initialize_transformer(nx, ny, &transformer)) {
    return 1;
  }
  if (0 != initialize_flow_field(&domain, &flow_field)) {
    return 1;
  }
  if (0 != initialize_flow_solver(&domain, &flow_solver)) {
    return 1;
  }
  if (0 != impose_initial_condition(&domain, transformer, &flow_field)) {
    return 1;
  }
  {
    double max_divergence = 0.;
    if (0 != check_divergence(&domain, &flow_field, &max_divergence)) {
      return 1;
    }
    printf("max divergence of initial flow field: % .1e\n", max_divergence);
  }
  for (;;) {
    const double time_max = 4e1;
    const double output_freq = 1e-1;
    static double time = 0.;
    static double output_next = output_freq;
    double dt = 0.;
    if (0 != integrate(&domain, transformer, &flow_solver, &flow_field, &dt)) {
      return 1;
    }
    time += dt;
    if (output_next < time) {
      static size_t output_counter = 0;
      char file_name[128] = {'\0'};
      snprintf(file_name, 128, OUTPUT_DIRECTORY "vorticity%010zu.npy", output_counter);
      output_vorticity(file_name, &domain, transformer, &flow_field);
      double max_divergence = 0.;
      if (0 != check_divergence(&domain, &flow_field, &max_divergence)) {
        return 1;
      }
      printf("time % .1e dt % .1e max divergence: % .1e\n", time, dt, max_divergence);
      output_counter += 1;
      output_next += output_freq;
    }
    if (time_max < time) {
      break;
    }
  }
  if (0 != finalize_domain(&domain)) {
    return 1;
  }
  if (0 != finalize_transformer(&transformer)) {
    return 1;
  }
  if (0 != finalize_flow_field(&flow_field)) {
    return 1;
  }
  if (0 != finalize_flow_solver(&flow_solver)) {
    return 1;
  }
  return 0;
}
