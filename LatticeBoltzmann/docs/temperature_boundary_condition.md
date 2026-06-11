# Derivation

We aim to compute

```math
g_2, g_5, g_6
```

at the bottom boundary of the domain using the remaining distribution function components.

We assume that the distribution function on the wall is at equilibrium, meaning it satisfies

```math
g_i
=
g_i^{eq}
\equiv
T^{\prime}
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
\right]
```

for $i = 2, 5, 6$.
Our goal is to determine the primed variables so that these relations are satisfied.

To this end, we start from the definitions of the temperature:

```math
T^w
=
g_0 + g_1 + g_2 + g_3 + g_4 + g_5 + g_6 + g_7 + g_8
```

and of the heat flux component:

```math
0
=
g_1 - g_3 + g_5 - g_6 - g_7 + g_8,
```

where we assume the wall is stationary.

From the equilibrium expression and the first relation, we get

```math
g_2 + g_5 + g_6
=
\left(
    w_2
    +
    w_5
    +
    w_6
\right)
T^{\prime}
=
T^w 
-
\left(
    g_0 + g_1 + g_3 + g_4 + g_7 + g_8
\right),
```

or equally

```math
T^{\prime}
=
\frac{
    T^w 
    -
    \left(
        g_0 + g_1 + g_3 + g_4 + g_7 + g_8
    \right)
}{
    W_1
    +
    2
    W_2
}.
```

Likewise, from the equilibrium expression and the second relation, we obtain

```math
g_5
-
g_6
=
\frac{2 W_2}{c_s^2}
u_x^{\prime}
T^{\prime}
=
-
\left(
    g_1
    -
    g_3
    -
    g_7
    +
    g_8
\right),
```

and thus

```math
u_x^{\prime}
=
-
\frac{
    c_s^2
}{
    2 W_2
}
\frac{
    g_1
    -
    g_3
    -
    g_7
    +
    g_8
}{
    T^{\prime}
}.
```

# Reference

- Inamuro et al., J. Comput. Phys., 2002

