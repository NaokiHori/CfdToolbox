#if !defined(ARRAY_H)
#define ARRAY_H

#include <stddef.h> // size_t

extern int prepare_array(
    const size_t nitems[NDIMS],
    double *** const array
);

extern int destroy_array(
    double *** const array
);

#endif // ARRAY_H
