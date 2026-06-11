#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "d2q9.h"
#include "domain.h"
#include "./fluid.h"
#include "./halo.h"
#include "./monitor.h"
#include "./output.h"
#include "./temperature.h"

typedef struct {
  double fluid;
  double temperature;
} diffusivities_t;

static diffusivities_t compute_diffusivities(
    const double ra,
    const double pr,
    const double ma
) {
  const double l = compute_domain_height();
  const double nu = (ma * sqrt(CS2)) * l * sqrt(pr / ra);
  const double kappa = nu / pr;
  diffusivities_t diffusivities = {
    .fluid = nu,
    .temperature = kappa,
  };
  return diffusivities;
}

static double compute_acceleration(
    const double ma
) {
  const double l = compute_domain_height();
  return pow(ma * sqrt(CS2), 2.) / l;
}

/**
 * time-step size in physical time units
 */
static double compute_time_step_size(
    const double acceleration
) {
  const double l = compute_domain_height();
  return sqrt(acceleration / l);
}

int main(
    void
) {
  // distribution functions
  df_t * const dfs[2] = {
    malloc(sizeof(df_t)),
    malloc(sizeof(df_t)),
  };
  if (NULL == dfs[0] || NULL == dfs[1]) {
    return 1;
  }
  // initialize distribution functions
  if (0 != initialize_fluid_distribution_function(dfs[0])) {
    return 1;
  }
  if (0 != initialize_temperature_distribution_function(dfs[1])) {
    return 1;
  }
  // prepare buffers to store macroscopic fields
  scalar_t * const density = malloc(sizeof(scalar_t));
  vector_t * const velocity = malloc(sizeof(vector_t));
  scalar_t * const temperature = malloc(sizeof(scalar_t));
  if (NULL == density) {
    return 1;
  }
  if (NULL == velocity) {
    return 1;
  }
  if (NULL == temperature) {
    return 1;
  }
  // parameters
  const double ra = 1e+8;
  const double pr = 1e+0;
  const double ma = 1e-1;
  // diffusivities and buoyancy acceleration
  const double l = compute_domain_height();
  const diffusivities_t diffusivities = compute_diffusivities(ra, pr, ma);
  const double acceleration = compute_acceleration(ma);
  // integrate in time
  // NOTE: in physical time units
  const double time_max = 5e+1;
  const double output_rate = 5e+0;
  const double monitor_rate = 5e-1;
  const double dt = compute_time_step_size(acceleration);
  for (;;) {
    static size_t step = 0;
    static double time = 0.;
    static double output_next = output_rate;
    static double monitor_next = monitor_rate;
    if (0 != compute_fluid_macroscopic_field(
        (const df_t *)dfs[0],
        density,
        velocity
    )) {
      return 1;
    }
    if (0 != compute_temperature_macroscopic_field(
        (const df_t *)dfs[1],
        temperature
    )) {
      return 1;
    }
    if (0 != process_fluid_collision(
        diffusivities.fluid,
        acceleration,
        (const scalar_t *)density,
        (const vector_t *)velocity,
        (const scalar_t *)temperature,
        dfs[0]
    )) {
      return 1;
    }
    if (0 != process_temperature_collision(
        diffusivities.temperature,
        (const vector_t *)velocity,
        (const scalar_t *)temperature,
        dfs[1]
    )) {
      return 1;
    }
    if (0 != process_fluid_streaming(dfs[0])) {
      return 1;
    }
    if (0 != process_temperature_streaming(dfs[1])) {
      return 1;
    }
    step += 1;
    time += dt;
    if (monitor_next < time) {
      double max_velocity = 0.;
      double nusselt = 0.;
      if (0 != monitor(
          l,
          diffusivities.temperature,
          (const vector_t *)velocity,
          (const scalar_t *)temperature,
          &max_velocity,
          &nusselt
      )) {
        return 0;
      }
      printf(
          "step %zu time % .1e Ma(max) % .1e Nu % .15e\n",
          step, time, max_velocity / sqrt(CS2), nusselt
      );
      monitor_next += monitor_rate;
    }
    if (output_next < time) {
#define FILE_NAME_MAX_LENGTH 128
      char file_name[FILE_NAME_MAX_LENGTH] = {'\0'};
      if (snprintf(file_name, FILE_NAME_MAX_LENGTH, "output/f%010zu.npy", step) < 0) {
        return 1;
      }
      if (0 != output_distribution_function(
          file_name,
          (const df_t *)dfs[0]
      )) {
        return 1;
      }
      if (snprintf(file_name, FILE_NAME_MAX_LENGTH, "output/g%010zu.npy", step) < 0) {
        return 1;
      }
      if (0 != output_distribution_function(
          file_name,
          (const df_t *)dfs[1]
      )) {
        return 1;
      }
      output_next += output_rate;
    }
    if (time_max < time) {
      break;
    }
  }
  free(dfs[0]);
  free(dfs[1]);
  free(density);
  free(velocity);
  free(temperature);
  return 0;
}

