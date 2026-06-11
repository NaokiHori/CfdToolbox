
.. include:: /reference.txt

##################
Governing Equation
##################

We consider a two-dimensional domain where a circular solid particle whose radius is :math:`\dim{a}` is positioned at the center.
Hereafter :math:`\dim{q}` denotes that the variable :math:`q` is dimensional (i.e., not non-dimensionalized), and vice versa.

The motion of liquid is governed by the incompressibility constraint:

.. math::

    \pder{}{\dim{u}_i}{\dim{x}_i}
    =
    0

as well as the momentum balance for Newtonian liquids:

.. math::

    \dim{\rho}
    \left(
        \pder{}{\dim{u}_i}{\dim{t}}
        +
        \dim{u}_j
        \pder{}{\dim{u}_i}{\dim{x}_j}
    \right)
    =
    \pder{}{\dim{\sigma}_{ij}}{\dim{x}_j}
    =
    -
    \pder{}{\dim{p}}{\dim{x}_i}
    +
    \dim{\mu}
    \pder{}{}{\dim{x}_j}
    \pder{}{\dim{u}_i}{\dim{x}_j}.

The particle emits a solute, which is transported into the liquid.
In particular, the concentration of the solute (mass per unit area) follows the advection-diffusion equation:

.. math::

    \pder{}{\dim{c}}{\dim{t}}
    +
    \dim{u}_j
    \pder{}{\dim{c}}{\dim{x}_j}
    =
    \dim{D}
    \pder{}{}{\dim{x}_j}
    \pder{}{\dim{c}}{\dim{x}_j},

where :math:`\dim{D}` is the diffusivity of the solute.

On the particle surface, we request the emission rate of the solute to be constant:

.. math::

    \dim{E}
    \equiv
    -
    \dim{D}
    \vat{
        \pder{}{\dim{c}}{\dim{\vr}}
    }{\dim{r} = \dim{a}}
    =
    const.,

which serves as a radial boundary condition with respect to the concentration field.
In the tangential direction, a slip velocity is prescribed:

.. math::

    \vat{
        \dim{\ut}
    }{\dim{r} = \dim{a}}
    =
    \dim{M}
    \vat{
        \frac{1}{\dim{\vr}}
        \pder{}{\dim{c}}{\dim{\vt}}
    }{\dim{r} = \dim{a}},

where :math:`\dim{M}` is called mobility.

The particle experiences hydrodynamic force exerted by the surround liquid, which is governed by the Newton-Euler equation:

.. math::

    &
    \dim{\rho} \pi \dim{a}^2 \der{}{\dim{v}_i}{\dim{t}}
    =
    \oint \dim{\sigma}_{ij} n_j d\dim{l},

    &
    \frac{1}{2} \dim{\rho} \pi \dim{a}^4 \der{}{\dim{\omega}_i}{\dim{t}}
    =
    \oint \epsilon_{ijk} \dim{r}_j \dim{\sigma}_{kl} n_l d\dim{l},

where we assume the particle is neutrally-buoyant.
:math:`\dim{v}_i` and :math:`\dim{\omega}_i` denote the translational and rotational velocity, and the right-hand-side integrals are applied to the interface along the particle.

Normalized equations are

.. math::

    \pder{}{u_i}{x_i}
    =
    0,

.. math::

    \rho
    \left(
        \pder{}{u_i}{t}
        +
        u_j
        \pder{}{u_i}{x_j}
    \right)
    =
    \frac{\dim{\mu} \dim{D}}{\dim{\rho} \dim{E} \dim{M} \dim{a}}
    \left(
        -
        \pder{}{p}{x_i}
        +
        \pder{}{}{x_j}
        \pder{}{u_i}{x_j}
    \right),

.. math::

    \pder{}{c}{t}
    +
    u_j
    \pder{}{c}{x_j}
    =
    \frac{\dim{\rho} \dim{D}}{\dim{E} \dim{a}}
    \frac{\dim{D}}{\dim{\rho} \dim{M}}
    \pder{}{}{x_j}
    \pder{}{c}{x_j},

with the corresponding boundary conditions on the particle:

.. math::

    \vat{
        \pder{}{c}{\vr}
    }{r = 1}
    =
    -
    \frac{\dim{E} \dim{a}}{\dim{\rho} \dim{D}},

.. math::

    \vat{
        \ut
    }{r = 1}
    =
    \frac{\dim{\rho} \dim{D}}{\dim{E} \dim{a}}
    \vat{
        \frac{1}{\vr}
        \pder{}{c}{\vt}
    }{r = 1}.

The non-dimensional Newton-Euler equation leads to

.. math::

    &
    \rho \pi a^2 \der{}{v_i}{t}
    =
    \frac{\dim{\mu} \dim{D}}{\dim{\rho} \dim{E} \dim{M} \dim{a}}
    \oint \sigma_{ij} n_j dl,

    &
    \frac{1}{2} \rho \pi a^4 \der{}{\omega_i}{t}
    =
    \frac{\dim{\mu} \dim{D}}{\dim{\rho} \dim{E} \dim{M} \dim{a}}
    \oint \epsilon_{ijk} r_j \sigma_{kl} n_l dl.

Note that we adopt :math:`\dim{a}` and :math:`\dim{E} \dim{M} / \dim{D}` as reference length and velocity scales to obtain the non-dimensional relations (|HU2019|).

In this project, we assume that the Reynolds number is sufficiently small so that the inertial effects on the momentum transport can be neglected.
By fixing the non-dimensional flux to be unity (which does not loss the generality), we obtain a set of equations governing the whole system:

.. math::

    \pder{}{u_i}{x_i}
    =
    0,

.. math::

    0_i
    =
    -
    \pder{}{p}{x_i}
    +
    \pder{}{}{x_j}
    \pder{}{u_i}{x_j},

.. math::

    \pder{}{c}{t}
    +
    u_j
    \pder{}{c}{x_j}
    =
    \frac{1}{Pe}
    \pder{}{}{x_j}
    \pder{}{c}{x_j},

with the corresponding boundary conditions on the particle:

.. math::

    \vat{
        \pder{}{c}{\vr}
    }{r = 1}
    =
    -
    1,

.. math::

    \vat{
        \ut
    }{r = 1}
    =
    \vat{
        \frac{1}{\vr}
        \pder{}{c}{\vt}
    }{r = 1}.

Note that the Newton-Euler equation in the limit of :math:`Re \rightarrow 0` leads to

.. math::

    0_i
    =
    \oint \sigma_{ij} n_j dl,

    0_i
    =
    \oint \epsilon_{ijk} r_j \sigma_{kl} n_l dl,

implying that the swimming particle (in the absence of external forces) is force-free and torque-free (|LAUGA2020|).

**************
Velocity Field
**************

Here we elaborate the velocity field.

To begin, we set a wall at :math:`\vr = R`.
The wall and the particle satisfy the impermeability condition:

.. math::

    \vat{\ur}{\vr = 1}
    =
    \vat{\ur}{\vr = R}
    =
    0.

On the inner wall, again, a slip velocity is enforced:

.. math::

    \vat{\ut}{\vr = 1}
    =
    \vat{
        \frac{1}{r}
        \pder{}{c}{\vt}
    }{\vr = 1},

whereas the outer wall follows the no-slip condition:

.. math::

    \vat{\ut}{\vr = R}
    =
    0.

Since the Stokes equation is a pure boundary-value problem, the entire flow field is uniquely determined once the boundary conditions are given.
In particular, the solution can be conveniently expressed using the stream function :math:`\psi \left( \vr, \vt, t \right)`, given by:

.. math::

    \psi \left( \vr, \vt, t \right)
    =
    \sum_k
    \Psi_k \left( \vr, t \right)
    \expp{I k \vt},

where :math:`I` is the imaginary unit, and :math:`\Psi_k` is the stream function in the frequency space.
Because of :math:`\psi \in \mathbb{R}`, :math:`\Psi_k` satisfies the complex-conjugate property, and thus it is sufficient to consider :math:`k \ge 0`.

As derived in :ref:`the appendix <appendix_stream_function>`, this is given by

.. math::

    \Psi_k
    =
    \begin{cases}
        k = 0
        &
        \text{unused}, \\
        k = 1
        &
        A_1 \vr^{-1}
        +
        B_1 \vr
        +
        E_1 \vr \log \vr
        +
        D_1 \vr^3, \\
        \text{otherwise}
        &
        A_k \vr^{-k}
        +
        B_k \vr^k
        +
        C_k \vr^{2 - k}
        +
        D_k \vr^{2 + k},
    \end{cases}

where the coefficients :math:`A_k, B_k, C_k, D_k, E_1` are determined using the boundary conditions at :math:`\vr = 1` and :math:`\vr = R`.

Assuming that :math:`R` is infinity, we request :math:`B_k, E_1, D_k` to be zero so that the solution converges, leading to

.. math::

    \Psi_k
    =
    A_k \vr^{-k}
    +
    C_k \vr^{2 - k}.

By imposing the boundary conditions at :math:`\vr = 1`:

.. math::

    \begin{pmatrix}
        1 & 1 \\
        - k & 2 - k \\
    \end{pmatrix}
    \begin{pmatrix}
        A_k \\
        C_k
    \end{pmatrix}
    =
    \begin{pmatrix}
        0 \\
        - \left( U_{\vt}^s \right)_k,
    \end{pmatrix}

we obtain

.. math::

    \Psi_k \left( \vr, t \right)
    =
    \frac{1 - \vr^2}{2 \vr^{k}}
    \left( U_{\vt}^s \right)_k.

In the above equation, :math:`U_{\vt}^s` is the azimuthal velocity along the particle surface, expressed in the frequency domain:

.. math::

    \left( U_{\vt}^s \right)_k
    =
    I k C_k^s,

or in the physical space:

.. math::

    \vat{\ut}{r = 1}
    =
    \vat{
        \frac{1}{\vr}
        \pder{}{\ut}{\vt}
    }{
        r = 1
    }.

:math:`C_k^s` represents the concentration field on the particle surface (at :math:`\vr = 1`) in the frequency domain:

.. math::

    c \left( \vr = 1, \vt, t \right)
    =
    \sum_k
    C_k^s \expp{I k \vt}.

*******************
Concentration Field
*******************

In addition to the constant-flux condition on the particle:

.. math::

    \vat{\pder{}{c}{\vr}}{\vr = 1}
    =
    -1,

we impose the Dirichlet condition on the outer boundary:

.. math::

    \vat{c}{\vr = R}
    =
    0.

In practice, we set :math:`R` to a finite value, which mismatches with the obtained stream function derived under the assumption that :math:`R \to \infty`.
This contradiction is originated by the fact that no analytical solution is available due to the non-linear term in the advection-diffusion equation.
We assume, nevertheless, that the approximation remains sufficiently accurate for large values of :math:`R`.

