#if !defined(FLOW_FIELD_H)
#define FLOW_FIELD_H

#include <complex.h>
#include "domain.h"
#include "transform.h"

typedef struct {
  double complex * ux;
  double complex * uy;
} flow_field_t;

extern int initialize_flow_field(
    const domain_t * const domain,
    flow_field_t * const flow_field
);

extern int finalize_flow_field(
    flow_field_t * const flow_field
);

extern int integrate_flow_field(
    const domain_t * const domain,
    transformer_t * const transformer,
    flow_field_t * const flow_field
);

#endif // FLOW_FIELD_H
