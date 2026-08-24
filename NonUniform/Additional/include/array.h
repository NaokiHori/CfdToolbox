#if !defined(ARRAY_H)
#define ARRAY_H

#include "domain.h"

typedef struct {
  double buffer[NY + 2][NX + 2];
} array_t;

#endif // ARRAY_H
