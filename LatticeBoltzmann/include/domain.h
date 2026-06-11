#if !defined(DOMAIN_H)
#define DOMAIN_H

#include "d2q9.h"

#define NX 256
#define NY 127

// data types to describe distribution functions and macroscopic fields
// NOTE: "+2" implies Q9
typedef double df_t[NY + 2][NX + 2][NDIRS];
typedef double scalar_t[NY + 2][NX + 2];
typedef double vector_t[NY + 2][NX + 2][NDIMS];

extern double compute_domain_height(
    void
);

#endif // DOMAIN_H
