#if !defined(UPDATE_CONCENTRATION_FIELD_H)
#define UPDATE_CONCENTRATION_FIELD_H

#include <stddef.h>
#include "domain.h"
#include "rdft.h"
#include "time_marcher.h"
#include "tri_diagonal_matrix.h"

extern int update_concentration_field(
    const double peclet,
    const domain_t * const domain,
    rdft_plan_t * const rdft_plan,
    tri_diagonal_matrix_t * const tri_diagonal_matrix,
    const double dt,
    const size_t runge_kutta_step,
    runge_kutta_buffer_t * const runge_kutta_buffer,
    double ** const concentration
);

#endif // UPDATE_CONCENTRATION_FIELD_H
