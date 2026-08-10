#include <math.h>
#include "./eigendecomposition.h"

static int compute_shift(
    const double * const main_diagonals,
    const double * const sub_diagonals,
    const size_t n_end,
    double * const shift
) {
  if (n_end < 2) {
    return 1;
  }
  const double diagonals[2] = {
    main_diagonals[n_end - 2],
    main_diagonals[n_end - 1],
  };
  const double sub_diagonal = sub_diagonals[n_end - 2];
  const double half_difference = 0.5 * (diagonals[0] - diagonals[1]);
  const double sign_half_difference = half_difference < 0. ? -1. : 1.;
  *shift =
    diagonals[1]
    -
    pow(sub_diagonal, 2.) / (
        half_difference
        +
        sign_half_difference * hypot(sub_diagonal, half_difference)
    );
  return 0;
}

static int chase_bulge(
    const size_t nitems,
    const size_t n_end,
    const double shift,
    double * const main_diagonals,
    double * const sub_diagonals,
    double * const eigenvectors_transposed
) {
  // rotation angle is determined by these quantities
  double values[2] = {
    main_diagonals[0] - shift,
    sub_diagonals[0],
  };
  for (size_t i = 0; i < n_end - 1; i++) {
    // Givens rotation applied to i-th and (i+1)-th rows
    // BEFORE
    //   (i    )-th row: [ . . s m s . . . ]
    //   (i + 1)-th row: [ . . b s m s . . ]
    // AFTER
    //   (i    )-th row: [ . . s m s b . . ]
    //   (i + 1)-th row: [ . . . s m s . . ]
    // - m: main-diagonal
    // - s: sub-diagonal
    // - b: bulge
    const double r = hypot(values[0], values[1]);
    if (0 < i) {
      // affected by eliminating the bulge
      sub_diagonals[i - 1] = r;
    }
    const double c = 0. == r ? 1. : + values[0] / r;
    const double s = 0. == r ? 0. : - values[1] / r;
    const double local_main_diagonals[2] = {
      main_diagonals[i    ],
      main_diagonals[i + 1],
    };
    const double local_sub_diagonal = sub_diagonals[i];
    // update the 2x2 diagonal block this rotation acts on
    // Q^T A Q = G A G^T =
    //   [  c -s ] [ d0  e ] [  c  s ]
    //   [  s  c ] [  e d1 ] [ -s  c ]
    const double cc = c * c;
    const double sc = s * c;
    const double ss = s * s;
    main_diagonals[i    ] = cc * local_main_diagonals[0] - 2. * sc * local_sub_diagonal + ss * local_main_diagonals[1];
    main_diagonals[i + 1] = ss * local_main_diagonals[0] + 2. * sc * local_sub_diagonal + cc * local_main_diagonals[1];
    sub_diagonals[i] = sc * local_main_diagonals[0] + (cc - ss) * local_sub_diagonal - sc * local_main_diagonals[1];
    // update (transposed) eigenvectors
    // V := V Q = V G^T
    // -> V^T := (V G^T)^T = G V^T
    for (size_t j = 0; j < nitems; j++) {
      // for each column (j), two rows (i and i + 1) are modulated
      double * const elements[2] = {
        eigenvectors_transposed + (i    ) * nitems + j,
        eigenvectors_transposed + (i + 1) * nitems + j,
      };
      const double values[2] = {
        *elements[0],
        *elements[1],
      };
      *elements[0] = + c * values[0] - s * values[1];
      *elements[1] = + s * values[0] + c * values[1];
    }
    if (n_end - 2 == i) {
      break;
    }
    // prepare for the next bulge chasing
    // [ b ] = [  c -s ] [ 0 ] = [ -s * e ]
    // [ s ] = [  s  c ] [ e ] = [  c * e ]
    values[0] = sub_diagonals[i];
    // new bulge
    values[1] = - s * sub_diagonals[i + 1];
    sub_diagonals[i + 1] *= c;
  }
  return 0;
}

static int build_identity_matrix(
    const size_t nitems,
    double * const matrix
) {
  for (size_t i = 0; i < nitems; i++) {
    for (size_t j = 0; j < nitems; j++) {
      matrix[i * nitems + j] = i == j ? 1. : 0.;
    }
  }
  return 0;
}

int compute_eigendecomposition(
    const size_t iter_max,
    const size_t nitems,
    double * const main_diagonals,
    double * const sub_diagonals,
    double * const eigenvalues,
    double * const eigenvectors_transposed
) {
  // qr algorithm
  build_identity_matrix(nitems, eigenvectors_transposed);
  // deflated position
  size_t n_end = nitems;
  for (size_t iter = 0; iter < iter_max; iter += 1) {
    double shift = 0.;
    compute_shift(main_diagonals, sub_diagonals, n_end, &shift);
    chase_bulge(nitems, n_end, shift, main_diagonals, sub_diagonals, eigenvectors_transposed);
    if (fabs(sub_diagonals[n_end - 2]) < 1e-15) {
      n_end -= 1;
    }
    if (n_end < 2) {
      for (size_t n = 0; n < nitems; n++) {
        eigenvalues[n] = main_diagonals[n];
      }
      return 0;
    }
  }
  // not converged
  return 1;
}

