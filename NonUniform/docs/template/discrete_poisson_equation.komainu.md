# Discrete Poisson equation

## Method

### 1D

For simplicity, let us first consider only the ${{ x }}$ direction. The one-dimensional discrete Poisson equation is

```math
{{
  index(lower_diagonal, i) * index(solution, i - 1)
  -
  (index(lower_diagonal, i) + index(upper_diagonal, i)) * index(solution, i)
  +
  index(upper_diagonal, i) * index(solution, i + 1)
  =
  index(right_hand_side, i)
}}.
```

We consider the eigenvalue decomposition of this equation.
Since it is a linear equation, it can be written as

```math
{{
  linear_operator * solution_vector
  =
  right_hand_side_vector
}}.
```

${{ linear_operator }}$ is a non-symmetric real tridiagonal matrix, which is to be decomposed into eigenvalues and eigenvectors.
Any real symmetric matrix ${{ real_symmetric_operator }}$ can be written with an orthogonal matrix ${{ orthogonal_operator }}$ (and its transpose ${{ transpose(orthogonal_operator) }}$) and a diagonal matrix ${{ diagonal_operator = diagonal(eigenvalue_vector) }}$ as

```math
{{
  real_symmetric_operator 
  =
  orthogonal_operator * diagonal_operator * transpose(orthogonal_operator)
}},
```

and we use this fact.
Since ${{ linear_operator }}$ is not symmetric, we need to convert it into a symmetric matrix in some way. From the structure of the original discrete Poisson equation, we find that it is enough to apply

```math
{{
  index_2d(symmetrization_operator, i, i)
  =
  sqrt(
    index((scale_factor), i)
  )
}},
```

that is,

```math
{{
  real_symmetric_operator
  =
  symmetrization_operator
  *
  linear_operator
  *
  inverse(symmetrization_operator)
}}.
```

Note that ${{ symmetrization_operator }}$ is a diagonal matrix, so its inverse ${{ inverse(symmetrization_operator) }}$ is obtained by simply taking the reciprocal of each element.
Putting these together, the target equation can be written as

```math
{{
  solution_vector
  =
  inverse(symmetrization_operator)
  *
  orthogonal_operator
  *
  inverse(diagonal_operator)
  *
  transpose(orthogonal_operator)
  *
  symmetrization_operator
  *
  right_hand_side_vector
}}.
```

In deriving this, we used the fact that ${{ orthogonal_operator }}$ is an orthogonal matrix, that is,

```math
{{ inverse(orthogonal_operator) = transpose(orthogonal_operator) }}.
```

Below, we summarize the procedure for actual numerical computation.

#### Step 1: Symmetrization

Compute ${{ real_symmetric_operator }}$.
The diagonal elements do not change. The elements just above and below the diagonal (which are equal because ${{ real_symmetric_operator }}$ is symmetric) are given by

```math
{{
  index_2d(real_symmetric_operator, i, i + 1)
}}
\equiv
{{
  index_2d(real_symmetric_operator, i + 1, i)
}}
=
{{
  fraction(
    1,
    sqrt(index((scale_factor), i))
  )
  *
  fraction(
    1,
    sqrt(index((scale_factor), i + 1))
  )
  *
  fraction(
    1,
    sqrt(index((scale_factor), i + fraction(1, 2)))
  )
}}.
```

#### Step 2: Eigenvalue decomposition

Compute the eigenvalue decomposition of the real symmetric (tridiagonal) matrix ${{ real_symmetric_operator }}$.
As a result, we obtain the eigenvalues ${{ diagonal_operator }}$ (a diagonal matrix) and the orthogonal matrix ${{ orthogonal_operator }}$ (whose columns are the eigenvectors).

#### Step 3: Compute the required matrices

Compute and store the following two square matrices of size ${{ dof }}$:

```math
{{ transpose(orthogonal_operator) * symmetrization_operator }},
```

```math
{{ inverse(symmetrization_operator) * orthogonal_operator }}.
```

Since ${{ symmetrization_operator }}$ is a diagonal matrix, this only requires multiplying or dividing the columns of ${{ transpose(orthogonal_operator) }}$ and the rows of ${{ orthogonal_operator }}$, respectively.

#### Step 4: Transform the right-hand side

By applying the transformation matrix, we move the right-hand side ${{ right_hand_side_vector }}$ into spectral space, ${{ spectral_right_hand_side_vector }}$:

```math
{{
  spectral_right_hand_side_vector
  =
  (transpose(orthogonal_operator) * symmetrization_operator) 
  *
  right_hand_side_vector
}}.
```

The matrix on the right-hand side was already computed in Step 3.

#### Step 5: Solve the Poisson equation

In spectral space, the Poisson equation is solved by simple division by the eigenvalues:

```math
{{
  spectral_solution_vector
  =
  inverse(diagonal_operator) * spectral_right_hand_side_vector
}},
```

that is,

```math
{{
  index(spectral_solution_vector, i)
  =
  fraction(
    index(spectral_right_hand_side_vector, i),
    index(eigenvalue_vector, i)
  )
}}.
```

#### Step 6: Return the solution to physical space

This is done by applying a transformation matrix (the inverse of the one used in Step 4):

```math
{{
  solution_vector
  =
  (inverse(symmetrization_operator) * orthogonal_operator) * spectral_solution_vector
}}.
```

Again, the matrix on the right-hand side was already computed in Step 3.

#### Summary

The key point is that Steps 1 to 3 only use the matrix on the left-hand side and do not use the right-hand side values at all.
The eigenvalue decomposition is a relatively heavy computation, but it only needs to be done once after the grid is generated. It does not need to be repeated while the fluid equations are being integrated.
When solving the Poisson equation inside the fluid solver, only Steps 4 to 6 are needed, and the operations performed at every time step are very simple.

### 2D

The method can be extended directly to two dimensions by applying a block matrix whose diagonal blocks are the matrix ${{ transpose(orthogonal_operator) }} {{ symmetrization_operator }}$.
The matrix notation may look like a big computation, but in a numerical code it is equivalent to performing a separate eigenvalue decomposition for each $j$.
Therefore, the extension to three dimensions is straightforward.

## Reference

- [Eigendecomposition of a matrix](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix)

