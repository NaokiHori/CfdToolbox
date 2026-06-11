#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constant.h"
#include "domain.h"
#include "memory.h"

// NOTE: fixed parameter, cannot be changed
static const double inner_radius = 1.;

int domain_initialize(
    const size_t nx,
    const size_t ny,
    const double outer_radius,
    domain_t * const domain
) {
  domain->nx = nx;
  domain->ny = ny;
  double ** const xf = &domain->xf;
  double ** const xc = &domain->xc;
  double ** const hxxf = &domain->hxxf;
  double ** const hxxc = &domain->hxxc;
  double ** const hyxf = &domain->hyxf;
  double ** const hyxc = &domain->hyxc;
  double ** const jdxf = &domain->jdxf;
  double ** const jdxc = &domain->jdxc;
  double * const dy = &domain->dy;
  *xf = memory_alloc(nx + 2, sizeof(double));
  *xc = memory_alloc(nx + 2, sizeof(double));
  *hxxf = memory_alloc(nx + 2, sizeof(double));
  *hxxc = memory_alloc(nx + 2, sizeof(double));
  *hyxf = memory_alloc(nx + 2, sizeof(double));
  *hyxc = memory_alloc(nx + 2, sizeof(double));
  *jdxf = memory_alloc(nx + 2, sizeof(double));
  *jdxc = memory_alloc(nx + 2, sizeof(double));
  // radial cell face / center positions
  // face positions are given by a quadratic function w.r.t. radial index,
  //   with the minimum grid size specified
  //   i = 1 and i = nx + 1 should agree with the inner and outer radii, respectively
  const double minimum_grid = 1. / nx;
  const double qa = 1. / nx * ((outer_radius - inner_radius) / nx - minimum_grid);
  const double qb = minimum_grid - 2. * qa;
  const double qc = inner_radius - qa - qb;
  (*xf)[0] = nan("");
  (*xf)[1] = inner_radius;
  (*xf)[nx + 1] = outer_radius;
  for (size_t i = 2; i <= nx; i++) {
    (*xf)[i] = qa * i * i + qb * i + qc;
  }
  (*xc)[0] = inner_radius;
  (*xc)[nx + 1] = outer_radius;
  for (size_t i = 1; i <= nx; i++) {
    (*xc)[i] = 0.5 * (*xf)[i] + 0.5 * (*xf)[i + 1];
  }
  // radial scale factors at radial cell face / center
  (*hxxf)[0] = nan("");
  for (size_t i = 1; i <= nx + 1; i++) {
    (*hxxf)[i] = (*xc)[i] - (*xc)[i - 1];
  }
  (*hxxc)[0] = nan("");
  (*hxxc)[nx + 1] = nan("");
  for (size_t i = 1; i <= nx; i++) {
    (*hxxc)[i] = (*xf)[i + 1] - (*xf)[i];
  }
  // azimuthal scale factors at radial cell face / center
  *dy = 2. * PI / ny;
  (*hyxf)[0] = nan("");
  for (size_t i = 1; i <= nx + 1; i++) {
    (*hyxf)[i] = (*xf)[i] * *dy;
  }
  (*hyxc)[0] = nan("");
  (*hyxc)[nx + 1] = nan("");
  for (size_t i = 1; i <= nx; i++) {
    (*hyxc)[i] = (*xc)[i] * *dy;
  }
  // jacobian determinants at radial cell face / center
  for (size_t i = 0; i <= nx + 1; i++) {
    (*jdxf)[i] = (*hxxf)[i] * (*hyxf)[i];
  }
  for (size_t i = 0; i <= nx + 1; i++) {
    (*jdxc)[i] = (*hxxc)[i] * (*hyxc)[i];
  }
  return 0;
}

int domain_finalize(
    domain_t * const domain
) {
  memory_free(domain->xf);
  memory_free(domain->xc);
  memory_free(domain->hxxf);
  memory_free(domain->hxxc);
  memory_free(domain->hyxf);
  memory_free(domain->hyxc);
  memory_free(domain->jdxf);
  memory_free(domain->jdxc);
  return 0;
}

