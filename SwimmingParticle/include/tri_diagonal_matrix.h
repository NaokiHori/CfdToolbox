#if !defined(TRI_DIAGONAL_MATRIX_H)
#define TRI_DIAGONAL_MATRIX_H

typedef struct {
  double ** l;
  double ** c;
  double ** u;
} tri_diagonal_matrix_t;

#endif // TRI_DIAGONAL_MATRIX_H
