#include "param.h"
#include "./compute_dux.h"
#include "./compute_dux/advx.h"
#include "./compute_dux/advy.h"
#include "./compute_dux/difx.h"
#include "./compute_dux/dify.h"
#include "./compute_dux/pres.h"

int compute_dux(
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const interface_field_t * const interface_field,
    const interface_solver_t * const interface_solver,
    const double dt,
    double ** const dux
) {
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      dux[j][i] = 0.;
    }
  }
  double ** const ux = flow_field->ux;
  double ** const uy = flow_field->uy;
  double ** const p = flow_field->p;
  const double c = 1. / Re;
  ux_advx(domain,     ux, dt, dux);
  ux_advy(domain, uy, ux, dt, dux);
  ux_difx(domain,  c, ux, dt, dux);
  ux_dify(domain,  c, ux, dt, dux);
  ux_pres(domain,      p, dt, dux);
  const double dx = domain->dx;
  double ** const vof = interface_field->vof;
  double ** const curvature = interface_solver->curvature;
#pragma omp parallel for
  for (size_t j = 1; j <= ny; j++) {
    for (size_t i = 2; i <= nx; i++) {
      dux[j][i] += dt / We
        * (0.5 * curvature[j][i - 1] + 0.5 * curvature[j][i])
        * (- vof[j][i - 1] + vof[j][i]) / dx;
    }
  }
  return 0;
}

