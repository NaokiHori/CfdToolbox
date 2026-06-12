#if !defined(INTERFACE_FIELD_H)
#define INTERFACE_FIELD_H

#include "domain.h" // domain_t

// diffusivity
extern const double interface_gamma;
// sharpness parameter
extern const double interface_epsilon;

typedef struct {
  double ** vof;
} interface_field_t;

extern int interface_field_init(
    const domain_t * const domain,
    interface_field_t * const interface_field
);

extern int interface_field_finalize(
    interface_field_t * const interface_field
);

#endif // INTERFACE_FIELD_H
