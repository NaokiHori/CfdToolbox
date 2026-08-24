#if !defined(FLOW_FIELD_H)
#define FLOW_FIELD_H

#include "array.h"
#include "domain.h" // domain_t

typedef struct {
  array_t * ux;
  array_t * uy;
  array_t * p;
  array_t * temperature;
  array_t * salinity;
} flow_field_t;

extern int flow_field_init(
    const domain_t * const domain,
    flow_field_t * const flow_field
);

extern int flow_field_finalize(
    flow_field_t * const flow_field
);

#endif // FLOW_FIELD_H
