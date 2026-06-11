# Derivation

We aim to compute

```math
f_2, f_5, f_6
```

at the bottom boundary of the domain using the remaining distribution function components.

We assume that the distribution function on the wall is at equilibrium, meaning it satisfies

```math
f_i
=
f_i^{eq}
\equiv
\rho^{\prime}
w_i
\left[
    1
    +
    \frac{
        \vec{c}_i
        \cdot
        \vec{u}^{\prime}
    }{
        c_s^2
    }
    +
    \frac{1}{2}
    \left(
        \frac{
            \vec{c}_i
            \cdot
            \vec{u}^{\prime}
        }{
            c_s^2
        }
    \right)^2
    -
    \frac{1}{2}
    \frac{
        \vec{u}^{\prime}
        \cdot
        \vec{u}^{\prime}
    }{
        c_s^2
    }
\right]
```

for $i = 2, 5, 6$.
Our goal is to determine the primed variables (the density and the stream-wise velocity component) so that these relations are satisfied.

To this end, we start from the definitions of the density:

```math
\rho^w
=
f_0 + f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8
```

and of the velocity components:

```math
\rho^w u_x^w
=
f_1 - f_3 + f_5 - f_6 - f_7 + f_8,
```

```math
\rho^w u_y^w
=
f_2 - f_4 + f_5 + f_6 - f_7 - f_8.
```

In this project, we assume that the walls are vertically stationary, giving $u_y^w = 0$, while $u_x^w$ is prescribed as the boundary condition.

By comparing the first relation with the condition $u_y^w = 0$, we find

```math
\rho^w
=
f_0 + f_1 + f_3 + 2 f_4 + 2 f_7 + 2 f_8,
```

where $\rho^w$, the density at the wall, is also an unknown.

From the equilibrium expression, we obtain

```math
\begin{aligned}
    f_2 + f_5 + f_6
    &
    =
    \rho^{\prime} 
    \left( w_2 + w_5 + w_6 \right) \\
    &
    =
    f_4 + f_7 + f_8
\end{aligned}
```

and

```math
\begin{aligned}
    f_5 - f_6
    &
    =
    \frac{
        \rho^{\prime}
        \left( w_5 + w_6 \right)
        u_x^{\prime}
    }{
        c_s^2
    } \\
    &
    =
    \rho^w u_x^w
    -
    f_1
    +
    f_3
    +
    f_7
    -
    f_8,
\end{aligned}
```

which yield

```math
\begin{aligned}
    \rho^{\prime} 
    &
    =
    \frac{
        f_4 + f_7 + f_8
    }{
        w_2 + w_5 + w_6
    } \\
    &
    =
    \frac{
        f_4 + f_7 + f_8
    }{
        W_1 + 2 W_2
    }
\end{aligned}
```

and

```math
\begin{aligned}
    u_x^{\prime}
    &
    =
    \frac{
        c_s^2
        \left(
            \rho^w u_x^w - f_1 + f_3 + f_7 - f_8
        \right)
    }{
        \rho^{\prime} \left( w_5 + w_6 \right)
    } \\
    &
    =
    \frac{
        c_s^2
        \left(
            \rho^w u_x^w - f_1 + f_3 + f_7 - f_8
        \right)
    }{
        2 \rho^{\prime} W_2
    }.
\end{aligned}
```

# Reference

- Inamuro et al., Phys. Fluids, 1995

