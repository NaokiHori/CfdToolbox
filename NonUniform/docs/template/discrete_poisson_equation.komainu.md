# Discrete Poisson equation

The two-dimensional Poisson equation discretized by the second-order accurate central-finite-difference schems in rectilinear coordinates leads to:

```math
{discrete_poisson_equation},
```

where

```math
{x_lower_diagonal_definition},
```

```math
{x_upper_diagonal_definition},
```

```math
{y_lower_diagonal_definition},
```

```math
{y_upper_diagonal_definition}.
```

This is a sparse linear system which looks like:

<p align="center">
  <img alt="sparse_matrix" src="./sparse_matrix.png" width="50%" />
</p>

where the red and blue cells contain non-zero and zero values, respectively.

To solve this system directly, we first aim to convert the original system to

```math
{discrete_poisson_equation_x_eigenspace},
```

which looks like:

<p align="center">
  <img alt="sparse_matrix_reduced" src="./sparse_matrix_reduced.png" width="50%" />
</p>

This conversion (eigendecomposition) is achieved by applying a block-diagonal operator ${block_diagonal_matrix_symbol}$:

```math
{block_diagonal_matrix_symbol}_{{ij}} {f}_{{jk}} {block_diagonal_matrix_symbol}_{{lk}}
=
{block_diagonal_matrix_symbol}_{{ij}} {g}_{{jk}} {block_diagonal_matrix_symbol}_{{lk}},
```

where ${block_diagonal_matrix_symbol}$ has ${orthogonal_matrix_symbol}$ (which is elaborated in the next part) repeatedly (for $N_y$ times) on the main diagonal.

The converted system is tridiagonal for each ${x_eigenvalue}$ and thus can be solved fairly easily by e.g., the Thomas algorithm.

## Spectral decomposition

To find the transform from ${f}$ to ${f_x_eigenspace}$, we focus on the one-dimensional equation:

```math
{discrete_poisson_equation_1d},
```

or equivalently:

```math
{discrete_poisson_equation_1d_matrix_form}.
```

Our objective is equivalent to decompose the linear operator ${laplace_operator}$.
To this end, we utilize that a real symmetric matrix can be decomposed as

```math
{real_symmetric_matrix_definition},
```

where ${orthogonal_matrix}$ and ${diagonal_matrix}$ are the orthogonal and diagonal matrices, respectively.

Note that ${laplace_operator}$ is not symmetric and thus the spectral decomposition is not directly applicable.
To symmetrize it, we use the following relation:

```math
{symmetrize_laplace_operator},
```

with ${symmetrization_matrix}$ being a diagonal matrix whose components are:

```math
{symmetrization_matrix_diagonal}.
```

As a consequence, the one-dimensional equation of interest leads to

```math
{discrete_poisson_equation_1d_matrix_form_decomposed},
```

or equivalently

```math
{discrete_poisson_equation_1d_matrix_form_decomposed_eigenspace},
```

where we use that inverse of an orthogonal matrix is equal to the transpose of it.

## Reference

- [Eigendecomposition of a matrix](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix)

