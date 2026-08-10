# Discrete Poisson equation

The two-dimensional Poisson equation discretized by the second-order accurate central-finite-difference schems in rectilinear coordinates leads to:

```math
{{\left( {l}_1 \right)_{i}} {{p}_{i-1, j}} - \left( {\left( {l}_1 \right)_{i}} + {\left( {u}_1 \right)_{i}} \right) {{p}_{i, j}} + {\left( {u}_1 \right)_{i}} {{p}_{i+1, j}} + {\left( {l}_2 \right)_{j}} {{p}_{i, j-1}} - \left( {\left( {l}_2 \right)_{j}} + {\left( {u}_2 \right)_{j}} \right) {{p}_{i, j}} + {\left( {u}_2 \right)_{j}} {{p}_{i, j+1}} = {{q}_{i, j}}},
```

where

```math
{{\left( {l}_1 \right)_{i}} = \frac{{1}}{{\left( {h}_{{\xi}^1} \right)_{i}}} \frac{{1}}{{\left( {h}_{{\xi}^1} \right)_{i-\frac{1}{2}}}}},
```

```math
{{\left( {u}_1 \right)_{i}} = \frac{{1}}{{\left( {h}_{{\xi}^1} \right)_{i}}} \frac{{1}}{{\left( {h}_{{\xi}^1} \right)_{i+\frac{1}{2}}}}},
```

```math
{{\left( {l}_2 \right)_{j}} = \frac{{1}}{{\left( {h}_{{\xi}^2} \right)_{j}}} \frac{{1}}{{\left( {h}_{{\xi}^2} \right)_{j-\frac{1}{2}}}}},
```

```math
{{\left( {u}_2 \right)_{j}} = \frac{{1}}{{\left( {h}_{{\xi}^2} \right)_{j}}} \frac{{1}}{{\left( {h}_{{\xi}^2} \right)_{j+\frac{1}{2}}}}}.
```

This is a sparse linear system which looks like:

![sparse_matrix](./sparse_matrix.png)

where the red and blue cells contain non-zero and zero values, respectively.

To solve this system directly, we first aim to convert the original system to

```math
{{\lambda_x} {{P}_{\lambda_k, j}} + {\left( {l}_2 \right)_{j}} {{P}_{\lambda_k, j-1}} - \left( {\left( {l}_2 \right)_{j}} + {\left( {u}_2 \right)_{j}} \right) {{P}_{\lambda_k, j}} + {\left( {u}_2 \right)_{j}} {{P}_{\lambda_k, j+1}} = {{Q}_{\lambda_k, j}}},
```

which is tri-diagonal for each ${\lambda_x}$ and thus can be solved fairly easily by e.g., the Thomas algorithm.

## Spectral decomposition

To find the transform from ${p}$ to ${P}$, we focus on the one-dimensional equation:

```math
{{\left( {l}_1 \right)_{i}} {{p}_{i-1}} - \left( {\left( {l}_1 \right)_{i}} + {\left( {u}_1 \right)_{i}} \right) {{p}_{i}} + {\left( {u}_1 \right)_{i}} {{p}_{i+1}} = {{q}_{i}}},
```

or equivalently:

```math
{{L_{ij}} {{p}_j} = {{q}_i}}.
```

Our objective is equivalent to decompose the linear operator ${L_{ij}}$.
To this end, we utilize that a real symmetric matrix can be decomposed as

```math
{{S_{ij}} = {Q_{ik}} {\Lambda_{kl}} {Q_{jl}}},
```

where ${Q_{ik}}$ and ${\Lambda_{kl}}$ are the orthogonal and diagonal matrices, respectively.

Note that ${L_{ij}}$ is not symmetric and thus the spectral decomposition is not directly applicable.
However, we use the following relation to diagonalize it:

```math
{{S_{ij}} = {D_{ik}} {L_{kl}} {D^{-1}_{lj}}},
```

with ${D_{ij}}$ being a diagonal matrix whose components are:

```math
{{D_{ii}} = {\sqrt{{\left( {h}_{{\xi}^1} \right)_{i}}}}}.
```

As a consequence, the one-dimensional equation of interest leads to

```math
{{\Lambda_{ik}} {Q_{lk}} {D_{lj}} {{p}_j} = {Q_{ki}} {D_{kj}} {{q}_j}},
```

or equivalently

```math
{{\Lambda_{ij}} {{P}_j} = {{Q}_i}},
```

where we use that inverse of an orthogonal matrix is equal to the transpose of it.

## Reference

- [Eigendecomposition of a matrix](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix)

