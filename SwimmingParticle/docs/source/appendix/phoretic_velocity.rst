
.. include:: /reference.txt

#################
Phoretic velocity
#################

According to |HU2019|, the phoretic velocity is given by

.. math::

    \vec{U}
    =
    -
    \vec{e}_{\vt}
    \frac{1}{2 \pi}
    \int_0^{2 \pi}
    M
    \vat{
        \frac{1}{\vr}
        \pder{}{c}{\vt}
    }{\vr = 1}
    d \vt
    =
    M
    \left[
        -
        \vec{e}_x
        \Re \left( C_1^s \right)
        +
        \vec{e}_y
        \Im \left( C_1^s \right)
    \right],

which is derived here.

Using

.. math::

    \vec{e}_{\vt}
    &
    =
    -
    \vec{e}_x \sin \vt
    +
    \vec{e}_y \cos \vt

    &
    =
    -
    \vec{e}_x \frac{I}{2} \left[
        \expp{I \vt}
        -
        \expp{- I \vt}
    \right]
    +
    \vec{e}_y \frac{1}{2} \left[
        \expp{I \vt}
        +
        \expp{- I \vt}
    \right]

and

.. math::

    \vat{\pder{}{c}{\vt}}{\vr = 1}
    =
    \sum_k I k C_k^s \expp{I k \vt},

we see that

.. math::

    -
    \frac{1}{2 \pi}
    \int_0^{2 \pi}
    M
    \frac{1}{\vr}
    \pder{}{c}{\vt}
    \vec{e}_{\vt}
    d \vt
    =
    &
    -
    \vec{e}_x
    \frac{1}{2 \pi}
    \int_0^{2 \pi}
    M
    \sum_k I k C_k^s \expp{I k \vt}
    \frac{I}{2} \left[
        \expp{I \vt}
        -
        \expp{- I \vt}
    \right]
    d \vt

    &
    -
    \vec{e}_y
    \frac{1}{2 \pi}
    \int_0^{2 \pi}
    M
    \sum_k I k C_k^s \expp{I k \vt}
    \frac{1}{2} \left[
        \expp{I \vt}
        +
        \expp{- I \vt}
    \right]
    d \vt

    =
    &
    +
    \vec{e}_x
    \frac{1}{2 \pi}
    \int_0^{2 \pi}
    M
    \sum_k
    \frac{k C_k^s}{2} \left[
        \expp{I \left( k + 1 \right) \vt}
        -
        \expp{I \left( k - 1 \right) \vt}
    \right]
    d \vt

    &
    -
    I
    \vec{e}_y
    \frac{1}{2 \pi}
    \int_0^{2 \pi}
    M
    \sum_k
    \frac{k C_k^s}{2} \left[
        \expp{I \left( k + 1 \right) \vt}
        +
        \expp{I \left( k - 1 \right) \vt}
    \right]
    d \vt.

It is readily apparent that the integral of the exponential functions in the range of :math:`\left[ 0, 2 \pi \right]` is zero.
Exceptions are when the exponent is zero:

.. math::

    \int_0^{2 \pi}
    1
    d \vt
    =
    2 \pi.

Thus the phoretic velocity leads to

.. math::

    \vec{U}
    =
    M
    \left[
        \vec{e}_x
        \left(
            -
            \frac{C_1^s}{2}
            -
            \frac{C_{-1}^s}{2}
        \right)
        +
        I
        \vec{e}_y
        \left(
            -
            \frac{C_1^s}{2}
            +
            \frac{C_{-1}^s}{2}
        \right)
    \right].

Because of :math:`c \in \mathbb{R}`, :math:`C_k^s` fulfills the complex-conjugate property:

.. math::

    \Re \left( C_{-1}^s \right)
    &
    =
    \Re \left( C_1^s \right),

    \Im \left( C_{-1}^s \right)
    &
    =
    -
    \Im \left( C_1^s \right).

Assigning these relations yields

.. math::

    \vec{U}
    =
    M
    \left[
        -
        \vec{e}_x
        \Re \left( C_1^s \right)
        +
        \vec{e}_y
        \Im \left( C_1^s \right)
    \right],

and thus we conclude

.. math::

    U_x
    =
    -
    M
    \Re \left( C_1^s \right),

    U_y
    =
    +
    M
    \Im \left( C_1^s \right).

