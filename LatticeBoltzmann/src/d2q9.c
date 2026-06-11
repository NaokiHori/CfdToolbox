#include "d2q9.h"

const double w[NDIRS] = {W0, W1, W1, W1, W1, W2, W2, W2, W2};
const int c[NDIRS][NDIMS] = {
  { 0,  0},
  {+1,  0},
  { 0, +1},
  {-1,  0},
  { 0, -1},
  {+1, +1},
  {-1, +1},
  {-1, -1},
  {+1, -1},
};

