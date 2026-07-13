#include <math.h>
#include <stdlib.h>
#include "array.h"
#include "param.h"
#include "particle_stuffs.h"

static int particle_init(
    const domain_t * const domain,
    particle_t * const particle
) {
  double * const density = &particle->density;
  double * const radius = &particle->radius;
  *density = 1.;
  *radius = 0.125 * fmin(domain->lx, domain->ly);
  particle->mass = *density * PI * pow(*radius, 2.);
  particle->moment_of_inertia = 0.5 * *density * PI * pow(*radius, 4.);
  particle->position.x = 0.25 * domain->lx;
  particle->position.y = 0.50 * domain->ly;
  particle->position.angle = nan(""); // angle is not tracked for simplicity
  particle->velocity.x = 0.;
  particle->velocity.y = 0.;
  particle->velocity.angle = 0.;
  return 0;
}

int particle_stuffs_init(
    const domain_t * const domain,
    particle_stuffs_t * const particle_stuffs
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  if (0 != particle_init(domain, &particle_stuffs->particle)) {
    return 1;
  }
  array_init(nx + 2, ny + 2, &particle_stuffs->acceleration_x);
  array_init(nx + 2, ny + 2, &particle_stuffs->acceleration_y);
  particle_stuffs->tolerance_threshold = 1e-12;
  return 0;
}

int particle_stuffs_finalize(
    particle_stuffs_t * const particle_stuffs
) {
  array_finalize(&particle_stuffs->acceleration_x);
  array_finalize(&particle_stuffs->acceleration_y);
  return 0;
}

