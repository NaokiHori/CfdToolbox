#if !defined(SAVE_H)
#define SAVE_H

#include <stddef.h> // size_t
#include "domain.h" // domain_t
#include "flow_field.h" // flow_field_t
#include "particle_stuffs.h" // particle_stuffs_t

extern int save(
    const size_t id,
    const size_t step,
    const double time,
    const domain_t * const domain,
    const flow_field_t * const flow_field,
    const particle_stuffs_t * const particle_stuffs
);

#endif // SAVE_H
