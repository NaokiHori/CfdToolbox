#include <stddef.h>
#include <stdlib.h>
#include "flow_field.h"

int initialize_flow_field(
    const domain_t * const domain,
    flow_field_t * const flow_field
) {
  const size_t knx = domain->knx;
  const size_t kny = domain->kny;
  flow_field->ux = malloc(knx * kny * sizeof(double complex));
  flow_field->uy = malloc(knx * kny * sizeof(double complex));
  return 0;
}

int finalize_flow_field(
    flow_field_t * const flow_field
) {
  free(flow_field->ux);
  free(flow_field->uy);
  return 0;
}

