#if !defined(LOAD_GRID_H)
#define LOAD_GRID_H

#define GRID_DIRECTORY "grid/"

#include <stdbool.h> // bool
#include <stddef.h> // size_t

extern int try_load_npy(
    const char file_name[],
    const size_t ndims,
    const size_t * const shape,
    const char dtype[],
    const bool is_fortran_order,
    const size_t element_size,
    void * const buffer
);

#endif // LOAD_GRID_H
