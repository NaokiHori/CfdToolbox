#include <stddef.h>
#include <stdio.h>
#include "array.h"
#include "domain.h"
#include "memory.h"
#include "particle.h"
#include "rdft.h"
#include "time_marcher.h"
#include "tri_diagonal_matrix.h"
#include "./compute_stream_function.h"
#include "./compute_surface_concentration.h"
#include "./compute_surface_velocity.h"
#include "./decide_time_step_size.h"
#include "./evaluate_right_hand_side.h"
#include "./initialize.h"
#include "./monitor.h"
#include "./output.h"
#include "./update_concentration_field.h"

typedef struct {
  const double freq;
  double next;
} schedule_t;

int main(
    void
) {
  const double peclet = 13.;
  domain_t domain = {0};
  domain_initialize(
      /* radial grids */    64,
      /* azimuthal grids */ 64,
      /* outer radius */    64.,
      &domain
  );
  // concentration field
  double ** concentration = NULL;
  prepare_array((size_t [NDIMS]){domain.nx + 2, domain.ny + 2}, &concentration);
  // stream function
  double ** stream_function = NULL;
  prepare_array((size_t [NDIMS]){domain.nx + 2, domain.ny + 2}, &stream_function);
  // particle velocity and position
  particle_t particle = {0};
  // intermediate buffers to use Runge-Kutta scheme
  runge_kutta_buffer_t runge_kutta_buffer = {0};
  prepare_array((size_t [NDIMS]){domain.nx + 2, domain.ny + 2}, &runge_kutta_buffer.non_linear);
  prepare_array((size_t [NDIMS]){domain.nx + 2, domain.ny + 2}, &runge_kutta_buffer.increment);
  // auxiliary buffer
  // - concentration on the object surface
  // - velocity on the object surface
  double * const surface_concentration = memory_alloc(domain.ny, sizeof(double));
  double * const surface_velocity = memory_alloc(domain.ny, sizeof(double));
  // auxiliary stuffs, FFTs and Thomas algorithm solver
  rdft_plan_t * rdft_plan = NULL;
  if (0 != rdft_init_plan(domain.ny, &rdft_plan)) {
    return 1;
  }
  tri_diagonal_matrix_t tri_diagonal_matrix = {0};
  prepare_array((size_t [NDIMS]){domain.ny + 2, domain.nx + 2}, &tri_diagonal_matrix.u);
  prepare_array((size_t [NDIMS]){domain.ny + 2, domain.nx + 2}, &tri_diagonal_matrix.c);
  prepare_array((size_t [NDIMS]){domain.ny + 2, domain.nx + 2}, &tri_diagonal_matrix.l);
  // impose initial condition
  initialize(&domain, concentration);
  particle.current_velocity[0] = 0.;
  particle.current_velocity[1] = 0.;
  particle.updated_velocity[0] = 0.;
  particle.updated_velocity[1] = 0.;
  particle.position[0] = 0.;
  particle.position[1] = 0.;
  // schedule events
  const double time_max = 1e+3 + 1e-1;
  schedule_t monitor_schedule = {
    .freq = 2.5e+0,
    .next = monitor_schedule.freq,
  };
  schedule_t output_schedule = {
    .freq = time_max / 100.,
    .next = output_schedule.freq,
  };
  // integrate in time
  for (;;) {
    static size_t step = 0;
    static double time = 0.;
    double dt = 0.;
    for (size_t runge_kutta_step = 0; runge_kutta_step < RUNGE_KUTTA_STEP_MAX; runge_kutta_step++) {
      // obtain velocity field (stream function)
      compute_surface_concentration(&domain, rdft_plan, concentration, surface_concentration);
      compute_surface_velocity(&domain, surface_concentration, surface_velocity);
      compute_stream_function(&domain, rdft_plan, surface_velocity, stream_function);
      // store current (old) translational velocity,
      //   which will be used to update position
      compute_particle_velocity(&domain, surface_concentration, particle.current_velocity);
      // time-step size is computed based on the computed velocity field
      if (0 == runge_kutta_step) {
        decide_time_step_size(peclet, &domain, stream_function, &dt);
      }
      // update concentration field
      evaluate_right_hand_side(peclet, &domain, runge_kutta_step, concentration, stream_function, &runge_kutta_buffer);
      update_concentration_field(peclet, &domain, rdft_plan, &tri_diagonal_matrix, dt, runge_kutta_step, &runge_kutta_buffer, concentration);
      // compute updated (new) translational velocity
      compute_surface_concentration(&domain, rdft_plan, concentration, surface_concentration);
      compute_particle_velocity(&domain, surface_concentration, particle.updated_velocity);
      // update position using old and new velocities
      update_particle_position(dt, runge_kutta_step, &particle);
    }
    step += 1;
    time += dt;
    if (monitor_schedule.next < time) {
      printf("step %8zu time % .2e time step size % .2e\n", step, time, dt);
      monitor(&domain, time, concentration, &particle);
      monitor_schedule.next += monitor_schedule.freq;
    }
    if (output_schedule.next < time) {
      output(&domain, concentration);
      output_schedule.next += output_schedule.freq;
    }
    if (time_max < time) {
      break;
    }
  }
  destroy_array(&concentration);
  destroy_array(&stream_function);
  destroy_array(&runge_kutta_buffer.non_linear);
  destroy_array(&runge_kutta_buffer.increment);
  destroy_array(&tri_diagonal_matrix.l);
  destroy_array(&tri_diagonal_matrix.c);
  destroy_array(&tri_diagonal_matrix.u);
  memory_free(surface_concentration);
  memory_free(surface_velocity);
  domain_finalize(&domain);
  if (0 != rdft_destroy_plan(&rdft_plan)) {
    return 1;
  }
  return 0;
}

