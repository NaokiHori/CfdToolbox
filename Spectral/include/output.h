#if !defined(OUTPUT_H)
#define OUTPUT_H

#include <stddef.h>

#define OUTPUT_DIRECTORY "output/"

extern int output(
    const char file_name[],
    const size_t ndims,
    const size_t * const nitems,
    const double * const buffer
);

#endif // OUTPUT_H
