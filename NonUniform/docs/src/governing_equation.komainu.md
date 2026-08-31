# Governing equation

The equations are expressed in the general orthogonal coordinate system:

```math
{incompressibility_constraint},
```

```math
{momentum_balance},
```

where ${scale_factor}$ is the scale factor:

```math
{scale_factor_definition}
```

and ${jacobian_determinant}$ is the Jacobian determinant:

```math
{jacobian_determinant_definition}.
```

Note that summation is taken with respect to the subscript $j$.
See [the other documentation](https://naokihori.github.io/SimpleTCSolver/index.html) for details.

The equations are integrated by means of the predictor-corrector approach: the velocity field is predicted:

```math
{predict_velocity},
```

which is corrected to enforce the incompressibility:

```math
{correct_velocity},
```

where ${scalar_potential}$ is scalar field obtained by solving the Poisson equation:

```math
{poisson_equation}.
```

Now, we write down this equation explicitly in a two-dimensional Cartesian coordinates.
The metrics have the following relations:

```math
{scale_factor_1}
=
{scale_factor_1} \left( {coordinate_1} \right),
```

```math
{scale_factor_2}
=
{scale_factor_2} \left( {coordinate_2} \right),
```

```math
{jacobian_determinant}
=
{scale_factor_1} {scale_factor_2}.
```

Since the scale factor in the first direction is independent in the second direction (and vice versa), the Poisson equation can be simplified as:

```math
{simplified_poisson_equation}.
```

