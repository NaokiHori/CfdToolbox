# Parameter configuration

The system is governed by two non-dimensional numbers: the Rayleigh number:

```math
Ra
\equiv
\frac{
    \beta g L^3 \left( \Delta T \right)
}{
    \nu \kappa
},
```

and the Prandtl number:

```math
Pr
\equiv
\frac{
    \nu
}{
    \kappa
}.
```

In the lattice Boltzmann method, one constraint is that the maximum velocity should be so small that it is much smaller than the sound speed $c_s \approx 0.58$ (in lattice units).
To satisfy this, the reference velocity given by:

```math
U_{ref}
=
\sqrt{
    \beta g L \left( \Delta T \right)
},
```

or equivalently the Mach number:

```math
Ma
\equiv
\frac{U_{ref}}{c_s}
```

should be small enough (to be prescribed).

With this, the gravitational acceleration and the kinematic viscosity are given by

```math
g
=
\frac{
    U_{ref}^2
}{
    \beta L \left( \Delta T \right)
},
```

and

```math
\nu
=
\frac{
    \sqrt{Pr}
}{
    \sqrt{Ra}
}
\sqrt{
    \beta g L^3 \left( \Delta T \right)
}
=
\frac{
    \sqrt{Pr}
}{
    \sqrt{Ra}
}
U_{ref}
L,
```

respectively.

Note that $L$ is determined by the number of grids in the wall-normal direction $N$:

```math
L = N + 1,
```

and we set the temperature difference between the walls $\Delta T$ to $1$ ($\beta \Delta T = 1$, strictly speaking) for convenience.

