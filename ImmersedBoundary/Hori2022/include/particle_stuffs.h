#if !defined(PARTICLE_STUFFS_H)
#define PARTICLE_STUFFS_H

#include "domain.h" // domain_t

// a general struct to store translational and rotational state
typedef struct {
  double x;
  double y;
  double angle;
} coordinate_t;


typedef struct {
  // metric
  double density; // assuming fluid density to be unity
  double radius;
  double mass;
  double moment_of_inertia;
  // gravity centered quantities
  coordinate_t position; // angle is not used
  coordinate_t velocity;
  coordinate_t position_increment; // angle is not used
  coordinate_t velocity_increment;
  coordinate_t integrated_acceleration;
  // fluid momentum inside particle computed at current / next steps
  coordinate_t internal_momentum_old;
  coordinate_t internal_momentum_new;
} particle_t;

typedef struct {
  particle_t particle;
  double ** acceleration_x;
  double ** acceleration_y;
  // determine particle convergence
  double tolerance_threshold;
} particle_stuffs_t;

extern int particle_stuffs_init(
    const domain_t * const domain,
    particle_stuffs_t * const particle_stuffs
);

extern int particle_stuffs_finalize(
    particle_stuffs_t * const particle_stuffs
);

#endif // PARTICLE_STUFFS_H
