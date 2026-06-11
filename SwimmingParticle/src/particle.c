#include "particle.h"
#include "time_marcher.h"

int compute_particle_velocity(
    const domain_t * const domain,
    double * const surface_concentration,
    double * velocity
) {
  const size_t ny = domain->ny;
  // the first mode gives the translational velocity
  //   of the object
  // normalize fft here
  velocity[0] = - surface_concentration[     1] / ny;
  velocity[1] = - surface_concentration[ny - 1] / ny;
  return 0;
}

int update_particle_position(
    const double dt,
    const size_t runge_kutta_step,
    particle_t * const particle
) {
  // update position following Crank-Nicolson scheme
  for (size_t dim = 0; dim < NDIMS; dim++) {
    const double a = RUNGE_KUTTA_COEFS.alpha[runge_kutta_step];
    const double old_v = particle->current_velocity[dim];
    const double new_v = particle->updated_velocity[dim];
    double * const p = &particle->position[dim];
    *p += a * dt * (0.5 * old_v + 0.5 * new_v);
  }
  return 0;
}

