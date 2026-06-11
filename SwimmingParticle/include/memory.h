#if !defined(MEMORY_H)
#define MEMORY_H

#include <stddef.h>

extern void * memory_alloc(
    const size_t nitems,
    const size_t size
);

extern void memory_free(
    void * const ptr
);

#endif // MEMORY_H
