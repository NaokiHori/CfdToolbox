#include <math.h>
#include "boundary_condition.h"
#include "exchange_halo.h"
#include "./immersed_boundary.h"

static double heaviside(
    const double x
) {
  return x <= 0. ? 0. : 1.;
}

static double compute_fraction(
    const double grid_size,
    const double radius,
    const double dx,
    const double dy
) {
  const double distance = fmax(
      sqrt(
          + pow(dx, 2.)
          + pow(dy, 2.)
      ),
      1e-16
  );
  const double nx = dx / distance;
  const double ny = dy / distance;
  const double lambda = fabs(nx) + fabs(ny);
  const double sigma = 0.05 * (1. - pow(lambda, 2.)) + 0.3;
  return 0.5 * (1. + tanh((radius - distance) / sigma / lambda / grid_size));
}

static size_t clip_i_of_ux(
    const size_t nx,
    const int i
) {
  if (i < 1) {
    return 1;
  }
  if (nx + 2 < (size_t)i) {
    return nx + 2;
  }
  return i;
}

static size_t clip_i_of_uy(
    const size_t nx,
    const int i
) {
  if (i < 1) {
    return 1;
  }
  if (nx + 1 < (size_t)i) {
    return nx + 1;
  }
  return i;
}

static size_t clip_j(
    const size_t n,
    const int j
) {
  if (j < 1) {
    return 1;
  }
  if (n + 1 < (size_t)j) {
    return n + 1;
  }
  return j;
}

static const int sub_min = 2;
static const int add_max = 3;

static int compute_i_extrema_of_ux(
    const size_t nx,
    const double dx,
    const double radius,
    const double x,
    size_t * const i_min,
    size_t * const i_max
) {
  // x = dx * (i - 1)
  // add / subtract enough to accommodate smoothed delta function
  const int index_min = floor((x - radius) / dx + 1.) - sub_min;
  const int index_max =  ceil((x + radius) / dx + 1.) + add_max;
  *i_min = clip_i_of_ux(nx, index_min);
  *i_max = clip_i_of_ux(nx, index_max);
  return 0;
}

static int compute_i_extrema_of_uy(
    const size_t nx,
    const double dx,
    const double radius,
    const double x,
    size_t * const i_min,
    size_t * const i_max
) {
  // x = dx * (i - 0.5)
  // add / subtract enough to accommodate smoothed delta function
  const int index_min = floor((x - radius) / dx + 0.5) - sub_min;
  const int index_max =  ceil((x + radius) / dx + 0.5) + add_max;
  *i_min = clip_i_of_uy(nx, index_min);
  *i_max = clip_i_of_uy(nx, index_max);
  return 0;
}

static int compute_j_extrema_of_ux(
    const size_t ny,
    const double dy,
    const double radius,
    const double y,
    size_t * const j_min,
    size_t * const j_max
) {
  // y = dy * (j - 0.5)
  // add / subtract enough to accommodate smoothed delta function
  const int index_min = floor((y - radius) / dy + 0.5) - sub_min;
  const int index_max =  ceil((y + radius) / dy + 0.5) + add_max;
  *j_min = clip_j(ny, index_min);
  *j_max = clip_j(ny, index_max);
  return 0;
}

static int compute_j_extrema_of_uy(
    const size_t ny,
    const double dy,
    const double radius,
    const double y,
    size_t * const j_min,
    size_t * const j_max
) {
  // y = dy * (j - 1)
  // add / subtract enough to accommodate smoothed delta function
  const int index_min = floor((y - radius) / dy + 1.) - sub_min;
  const int index_max =  ceil((y + radius) / dy + 1.) + add_max;
  *j_min = clip_j(ny, index_min);
  *j_max = clip_j(ny, index_max);
  return 0;
}

static double compute_acceleration(
    const double fluid_velocity,
    const double particle_velocity,
    const double dt,
    const double fraction
) {
  return heaviside(1. - fraction) * fraction * (particle_velocity - fluid_velocity) / dt;
}

int immersed_boundary_reset_particle_increments(
    particle_stuffs_t * const particle_stuffs
) {
  particle_t * const particle = &particle_stuffs->particle;
  coordinate_t * const position_increment = &particle->position_increment;
  coordinate_t * const velocity_increment = &particle->velocity_increment;
  position_increment->x = 0.;
  position_increment->y = 0.;
  position_increment->angle = nan(""); // not in use
  velocity_increment->x = 0.;
  velocity_increment->y = 0.;
  velocity_increment->angle = 0.;
  return 0;
}

int immersed_boundary_compute_internal_momentum(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    particle_stuffs_t * const particle_stuffs,
    const bool update_new
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double ly = domain->ly;
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double grid_size = fmin(dx, dy);
  const double cell_size = dx * dy;
  particle_t * const particle = &particle_stuffs->particle;
  const double radius = particle->radius;
  double ** const ux = flow_field->ux;
  double ** const uy = flow_field->uy;
  double internal_momentum_x = 0.;
  double internal_momentum_y = 0.;
  double internal_momentum_angle = 0.;
  for (int y_periodic_state = -1; y_periodic_state <= 1; y_periodic_state++) {
    const double center_x = particle->position.x + particle->position_increment.x;
    const double center_y = particle->position.y + particle->position_increment.y + ly * y_periodic_state;
    size_t i_min = 0;
    size_t i_max = 0;
    size_t j_min = 0;
    size_t j_max = 0;
    compute_i_extrema_of_ux(nx, dx, radius, center_x, &i_min, &i_max);
    compute_j_extrema_of_ux(ny, dy, radius, center_y, &j_min, &j_max);
#pragma omp parallel for reduction(+:internal_momentum_x, internal_momentum_angle)
    for (size_t j = j_min; j < j_max; j++) {
      const double y = dy * (j - 0.5);
      for (size_t i = i_min; i < i_max; i++) {
        const double x = dx * (i - 1.);
        const double rx = x - center_x;
        const double ry = y - center_y;
        const double fraction = compute_fraction(
            grid_size,
            radius,
            rx,
            ry
        );
        const double velocity = ux[j][i];
        internal_momentum_x += fraction * velocity * cell_size;
        internal_momentum_angle -= fraction * ry * velocity * cell_size;
      }
    }
    compute_i_extrema_of_uy(nx, dx, radius, center_x, &i_min, &i_max);
    compute_j_extrema_of_uy(ny, dy, radius, center_y, &j_min, &j_max);
#pragma omp parallel for reduction(+:internal_momentum_y, internal_momentum_angle)
    for (size_t j = j_min; j < j_max; j++) {
      const double y = dy * (j - 1.);
      for (size_t i = i_min; i < i_max; i++) {
        const double x = dx * (i - 0.5);
        const double rx = x - center_x;
        const double ry = y - center_y;
        const double fraction = compute_fraction(
            grid_size,
            radius,
            rx,
            ry
        );
        const double velocity = uy[j][i];
        internal_momentum_y += fraction * velocity * cell_size;
        internal_momentum_angle += fraction * rx * velocity * cell_size;
      }
    }
  }
  coordinate_t * const internal_momentum = update_new
    ? &particle->internal_momentum_new
    : &particle->internal_momentum_old;
  internal_momentum->x = internal_momentum_x;
  internal_momentum->y = internal_momentum_y;
  internal_momentum->angle = internal_momentum_angle;
  return 0;
}

int immersed_boundary_exchange_momentum(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    particle_stuffs_t * const particle_stuffs,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double ly = domain->ly;
  const double dx = domain->dx;
  const double dy = domain->dy;
  const double grid_size = fmin(dx, dy);
  const double cell_size = dx * dy;
  double ** const fluid_ux = flow_field->ux;
  double ** const fluid_uy = flow_field->uy;
  particle_t * const particle = &particle_stuffs->particle;
  const double radius = particle->radius;
  const double center_ux = particle->velocity.x + particle->velocity_increment.x;
  const double center_uy = particle->velocity.y + particle->velocity_increment.y;
  const double omega = particle->velocity.angle + particle->velocity_increment.angle;
  double ** const acceleration_x = particle_stuffs->acceleration_x;
  double ** const acceleration_y = particle_stuffs->acceleration_y;
  double particle_integrated_acceleration_x = 0.;
  double particle_integrated_acceleration_y = 0.;
  double particle_integrated_acceleration_angle = 0.;
  for (int y_periodic_state = -1; y_periodic_state <= 1; y_periodic_state++) {
    const double center_x = particle->position.x + particle->position_increment.x;
    const double center_y = particle->position.y + particle->position_increment.y + ly * y_periodic_state;
    size_t i_min = 0;
    size_t i_max = 0;
    size_t j_min = 0;
    size_t j_max = 0;
    compute_i_extrema_of_ux(nx, dx, radius, center_x, &i_min, &i_max);
    compute_j_extrema_of_ux(ny, dy, radius, center_y, &j_min, &j_max);
#pragma omp parallel for reduction(+:particle_integrated_acceleration_x, particle_integrated_acceleration_angle)
    for (size_t j = j_min; j < j_max; j++) {
      const double y = dy * (j - 0.5);
      for (size_t i = i_min; i < i_max; i++) {
        const double x = dx * (i - 1.);
        const double rx = x - center_x;
        const double ry = y - center_y;
        const double fraction = compute_fraction(grid_size, radius, rx, ry);
        const double particle_ux = center_ux - omega * ry;
        const double acceleration = compute_acceleration(fluid_ux[j][i], particle_ux, dt, fraction);
        particle_integrated_acceleration_x -= acceleration * cell_size;
        particle_integrated_acceleration_angle += ry * acceleration * cell_size;
        acceleration_x[j][i] = acceleration;
      }
    }
    compute_i_extrema_of_uy(nx, dx, radius, center_x, &i_min, &i_max);
    compute_j_extrema_of_uy(ny, dy, radius, center_y, &j_min, &j_max);
#pragma omp parallel for reduction(+:particle_integrated_acceleration_y, particle_integrated_acceleration_angle)
    for (size_t j = j_min; j < j_max; j++) {
      const double y = dy * (j - 1.);
      for (size_t i = i_min; i < i_max; i++) {
        const double x = dx * (i - 0.5);
        const double rx = x - center_x;
        const double ry = y - center_y;
        const double fraction = compute_fraction(grid_size, radius, rx, ry);
        const double particle_uy = center_uy + omega * rx;
        const double acceleration = compute_acceleration(fluid_uy[j][i], particle_uy, dt, fraction);
        particle_integrated_acceleration_y -= acceleration * cell_size;
        particle_integrated_acceleration_angle -= rx * acceleration * cell_size;
        acceleration_y[j][i] = acceleration;
      }
    }
  }
  particle->integrated_acceleration.x = particle_integrated_acceleration_x;
  particle->integrated_acceleration.y = particle_integrated_acceleration_y;
  particle->integrated_acceleration.angle = particle_integrated_acceleration_angle;
  return 0;
}

int immersed_boundary_correct_velocity(
    const domain_t * const domain,
    flow_field_t * const flow_field,
    particle_stuffs_t * const particle_stuffs,
    const double dt
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  double ** const ux = flow_field->ux;
  double ** const uy = flow_field->uy;
  double ** const acceleration_x = particle_stuffs->acceleration_x;
  double ** const acceleration_y = particle_stuffs->acceleration_y;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      ux[j][i] += acceleration_x[j][i] * dt;
    }
  }
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 1; i <= nx; i++) {
      uy[j][i] += acceleration_y[j][i] * dt;
    }
  }
  if (0 != impose_boundary_condition_ux_x(domain, ux)) {
    return 1;
  }
  if (0 != exchange_halo_y(domain, ux)) {
    return 1;
  }
  if (0 != impose_boundary_condition_uy_x(domain, uy)) {
    return 1;
  }
  if (0 != exchange_halo_y(domain, uy)) {
    return 1;
  }
  return 0;
}

int immersed_boundary_update_particle_increments(
    particle_stuffs_t * const particle_stuffs,
    const double dt,
    double * const tolerance
) {
  particle_t * const particle = &particle_stuffs->particle;
  const double mass = particle->mass;
  const double moment_of_inertia = particle->moment_of_inertia;
  // NOTE: two position increments are defined for clarity
  const coordinate_t * const integrated_acceleration = &particle->integrated_acceleration;
  const coordinate_t * const velocity_old = &particle->velocity;
  const coordinate_t * const position_increment_old = &particle->position_increment;
  const coordinate_t * const internal_momentum_old = &particle->internal_momentum_old;
  const coordinate_t * const internal_momentum_new = &particle->internal_momentum_new;
  coordinate_t * const position_increment = &particle->position_increment;
  coordinate_t * const velocity_increment = &particle->velocity_increment;
  // euler backward
  velocity_increment->x =
    + dt / mass * integrated_acceleration->x
    + 1. / mass * (internal_momentum_new->x - internal_momentum_old->x);
  velocity_increment->y =
    + dt / mass * integrated_acceleration->y
    + 1. / mass * (internal_momentum_new->y - internal_momentum_old->y);
  velocity_increment->angle =
    + dt / moment_of_inertia * integrated_acceleration->angle
    + 1. / moment_of_inertia * (internal_momentum_new->angle - internal_momentum_old->angle);
  const coordinate_t velocity_new = {
    .x = velocity_old->x + velocity_increment->x,
    .y = velocity_old->y + velocity_increment->y,
    .angle = nan(""), // not in use
  };
  // crank nicolson
  const coordinate_t position_increment_new = {
    .x = 0.5 * dt * (velocity_old->x + velocity_new.x),
    .y = 0.5 * dt * (velocity_old->y + velocity_new.y),
    .angle = nan(""), // not in use
  };
  *tolerance = sqrt(
      + pow(position_increment_new.x - position_increment_old->x, 2.)
      + pow(position_increment_new.y - position_increment_old->y, 2.)
  );
  *position_increment = position_increment_new;
  return 0;
}

int immersed_boundary_update_particle(
    const domain_t * const domain,
    particle_stuffs_t * const particle_stuffs
) {
  particle_t * const particle = &particle_stuffs->particle;
  double * const x = &particle->position.x;
  double * const y = &particle->position.y;
  particle->velocity.x += particle->velocity_increment.x;
  particle->velocity.y += particle->velocity_increment.y;
  particle->velocity.angle += particle->velocity_increment.angle;
  *x += particle->position_increment.x;
  *y += particle->position_increment.y;
  // correct periodicity
  const double ly = domain->ly;
  if (*y < 0.) {
    *y += ly;
  } else if (ly < *y) {
    *y -= ly;
  }
  return 0;
}

