#include <errno.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include "domain.h"
#include "./compute_surface_concentration.h"
#include "./monitor.h"

static int monitor_particle(
    const double time,
    const particle_t * const particle,
    const char fopen_mode[]
) {
  const char file_name[] = "output/particle.dat";
  errno = 0;
  FILE * const fp = fopen(file_name, fopen_mode);
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  fprintf(fp, "% .7e % .7e % .7e % .7e % .7e\n", time, particle->updated_velocity[0], particle->updated_velocity[1], particle->position[0], particle->position[1]);
  fclose(fp);
  return 0;
}

static int monitor_concentration_field(
    const domain_t * const domain,
    const double time,
    double ** const concentration,
    const char fopen_mode[]
) {
  const char file_name[] = "output/concentration.dat";
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  const double * const jdxc = domain->jdxc;
  double min = DBL_MAX;
  double max = - DBL_MAX;
  double sum = 0.;
  for (size_t i = 1; i <= nx; i++) {
    for (size_t j = 1; j <= ny; j++) {
      min = fmin(min, concentration[i][j]);
      max = fmax(max, concentration[i][j]);
      sum += concentration[i][j] * jdxc[i];
    }
  }
  errno = 0;
  FILE * const fp = fopen(file_name, fopen_mode);
  if (NULL == fp) {
    perror(file_name);
    return 1;
  }
  fprintf(fp, "% .7e % .7e % .7e % .7e\n", time, min, max, sum);
  fclose(fp);
  return 0;
}

int monitor(
    const domain_t * const domain,
    const double time,
    double ** const concentration,
    const particle_t * const particle
) {
  static size_t counter = 0;
  const char * const fopen_mode = 0 == counter ? "w" : "a";
  monitor_particle(time, particle, fopen_mode);
  monitor_concentration_field(domain, time, concentration, fopen_mode);
  counter += 1;
  return 0;
}

