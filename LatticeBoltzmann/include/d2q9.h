#if !defined(D2Q9_H)
#define D2Q9_H

#define NDIMS 2
#define NDIRS 9

#define W0 16. / 36.
#define W1  4. / 36.
#define W2  1. / 36.

extern const double w[NDIRS];
extern const int c[NDIRS][NDIMS];

// squared speed of sound (in lattice unit)
#define CS2 0.3333333333333333

#endif // D2Q9_H
