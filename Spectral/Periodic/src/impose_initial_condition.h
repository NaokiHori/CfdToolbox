#if !defined(IMPOSE_INITIAL_CONDITION_H)
#define IMPOSE_INITIAL_CONDITION_H

#include "domain.h"
#include "flow_field.h"
#include "transform.h"

extern int impose_initial_condition(
    const domain_t * const domain,
    transformer_t * const transformer,
    flow_field_t * const flow_field
);

#endif // IMPOSE_INITIAL_CONDITION_H
