# Governing equation

The equations are expressed in the general orthogonal coordinate system:

```math
{{
  1
  /
  jacobian_determinant
  *
  partial_operator(
    jacobian_determinant
    /
    scale_factor(j)
    *
    velocity(j),
    coordinate(j)
  )
  =
  0
}},
```

```math
{{
  partial(
    velocity(i),
    time
  )
  =
  -
  fraction(
    velocity(j),
    scale_factor(j)
  )
  *
  partial(
    velocity(i),
    coordinate(j)
  )
  -
  1
  /
  scale_factor(i)
  *
  partial(
    pressure,
    coordinate(i)
  )
  +
  fraction(
    1,
    reynolds
  )
  *
  fraction(
    1,
    jacobian_determinant
  )
  *
  partial_operator(
    fraction(
      jacobian_determinant,
      scale_factor(j)
    )
    *
    fraction(
      1,
      scale_factor(j)
    )
    *
    partial(
      velocity(i),
      coordinate(j)
    ),
    coordinate(j)
  )
}},
```

where ${{ scale_factor(i) }}$ is the scale factor:

```math
{{
  partial(
    cartesian(i),
    coordinate(i)
  )
}},
```

and ${{ jacobian_determinant }}$ is the Jacobian determinant defined as:

```math
{{
  product_over(i, scale_factor(i))
}}.
```

Note that summation is taken with respect to ${{ j }}$.
See [the other documentation](https://naokihori.github.io/SimpleTCSolver/index.html) for details.

The equations are integrated by means of the predictor-corrector approach: the velocity field is predicted:

```math
{{
  step_at(velocity(i), predicted)
  =
  step_at(velocity(i), current)
  +
  definite_integral(
    step_at(time, current),
    step_at(time, next),
    vector(right_hand_side, i),
    time
  )
}},
```

which is corrected to enforce the incompressibility:

```math
{{
  step_at(velocity(i), next)
  =
  step_at(velocity(i), predicted)
  -
  definite_integral(
    step_at(time, current),
    step_at(time, next),
    fraction(
      1,
      scale_factor(i)
    )
    *
    partial(
      scalar_potential,
      coordinate(i)
    ),
    time
  )
}},
```

where ${{ scalar_potential }}$ is scalar field obtained by solving the Poisson equation:

```math
{{
  definite_integral(
    step_at(time, current),
    step_at(time, next),
    fraction(
      1,
      jacobian_determinant
    )
    *
    partial_operator(
      fraction(
        jacobian_determinant,
        scale_factor(j)
      )
      *
      fraction(
        1,
        scale_factor(j)
      )
      *
      partial(
        scalar_potential,
        coordinate(j)
      ),
      coordinate(j)
    ),
    time
  )
  =
  fraction(
    1,
    jacobian_determinant
  )
  *
  partial_operator(
    fraction(
      jacobian_determinant,
      scale_factor(j)
    )
    *
    step_at(velocity(j), predicted),
    coordinate(j)
  )
}}.
```

Now, we consider this equation explicitly in a two-dimensional Cartesian coordinates.
Since each scale factor ${{ scale_factor(i) }}$ is a function of the coordinate in the corresponding direction ${{ cartesian(i) }}$, or in other words since the scale factor in the first direction is independent in the second direction (and vice versa), the Poisson equation can be simplified as:

```math
{{
  definite_integral(
    step_at(time, current),
    step_at(time, next),
    fraction(
      1,
      scale_factor(1)
    )
    *
    partial_operator(
      fraction(
        1,
        scale_factor(1)
      )
      *
      partial(
        scalar_potential,
        coordinate(1)
      ),
      coordinate(1)
    )
    +
    fraction(
      1,
      scale_factor(2)
    )
    *
    partial_operator(
      fraction(
        1,
        scale_factor(2)
      )
      *
      partial(
        scalar_potential,
        coordinate(2)
      ),
      coordinate(2)
    ),
    time
  )
  =
  fraction(
    1,
    scale_factor(1)
  )
  *
  partial(
    step_at(velocity(1), predicted),
    coordinate(1)
  )
  +
  fraction(
    1,
    scale_factor(2)
  )
  *
  partial(
    step_at(velocity(2), predicted),
    coordinate(2)
  )
}}.
```

