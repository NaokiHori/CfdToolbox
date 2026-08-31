#if !defined(BOUNDARY_CONDITION_H)
#define BOUNDARY_CONDITION_H

#include "array.h" // array_t

extern int impose_boundary_condition_ux_x(
    array_t * const ux_array
);

extern int impose_boundary_condition_ux_y(
    array_t * const ux_array
);

extern int impose_boundary_condition_uy_x(
    array_t * const uy_array
);

extern int impose_boundary_condition_uy_y(
    array_t * const uy_array
);

extern int impose_boundary_condition_scalar_x(
    array_t * const scalar_array
);

extern int impose_boundary_condition_scalar_y(
    array_t * const scalar_array
);

#endif // BOUNDARY_CONDITION_H
