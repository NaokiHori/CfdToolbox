We consider a solid circular object, whose position, translational velocity, and angular velocity of the gravity center are denoted by ${{\Xi}_i}$, ${{\Upsilon}_i}$, and ${{\Omega}_i}$, respectively.
Due to the rigid body motion, the velocity vector at any point of the object ${{X}_i}$ is described as:

```math
{{U}_i}
=
{{\Upsilon}_i}
+
{\epsilon_{ijk}} {{\Omega}_j} \left( {{X}_k} - {{\Xi}_k} \right).
```

The fluid velocity field is updated using the Euler-forward scheme, without incorporating the effect of the immersed object:

```math
{{{u}_i}^{*}}
=
{{{u}_i}^{n}}
+
{\Delta {t}}
\left(
    - {{{u}_j}^{n}} \frac{\partial {{{u}_i}^{n}}}{\partial {{x}_j}} + \frac{{1}}{{\rho^f}} \frac{\partial {{\sigma_{ij}}^{n}}}{\partial {{x}_j}} + {{g}_i}
\right).
```

To take into account the effect of the particle by enforcing the desired fluid velocity on the particle surface, fluid and particle momenta are exchanged through:

```math
{{{a}_i}^{IBM}}
=
{H \left( {1} - {{\phi}^{n}} \right)} {{\phi}^{n}} \frac{{{{U}_i}^{n}} - {{{u}_i}^{*}}}{{\Delta {t}}},
```

where ${\phi}$ is the local volume fraction of the particle.

The fluid velocity is further corrected to account for this effect:

```math
{{{u}_i}^{**}}
=
{{{u}_i}^{*}}
+
{{{a}_i}^{IBM}}
{\Delta {t}}.
```

The same variable ${{{a}_i}^{IBM}}$ is used to evaluate the particle acceleration:

```math
{{{{\Upsilon}_i}^{n + 1}} = {{{\Upsilon}_i}^{n}} - \frac{{1}}{{m^p}} \left( {\Delta {t}} {\int_{}^{} {\rho^f} {{{a}_i}^{IBM}} d{V}} + {\int_{}^{} {\rho^f} {{\phi}^{n + 1}} {{{u}_i}^{n + 1}} d{V}} - {\int_{}^{} {\rho^f} {{\phi}^{n}} {{{u}_i}^{n}} d{V}} \right)},
```

```math
{{{{\Omega}_i}^{n + 1}} = {{{\Omega}_i}^{n}} - \frac{{1}}{{I^p}} \left( {\Delta {t}} {\int_{}^{} {\rho^f} {\epsilon_{ijk}} \left( {{X}_j} - {{\Xi}_j} \right) {{{a}_k}^{IBM}} d{V}} + {\int_{}^{} {\rho^f} {{\phi}^{n + 1}} {\epsilon_{ijk}} \left( {{X}_j} - {{\Xi}_j} \right) {{{u}_k}^{n + 1}} d{V}} - {\int_{}^{} {\rho^f} {{\phi}^{n + 1}} {\epsilon_{ijk}} \left( {{X}_j} - {{\Xi}_j} \right) {{{u}_k}^{n}} d{V}} \right)},
```

where ${m^p}$ and ${I^p}$ are the mass and the moment of inertia of the particle:

```math
{m^p}
=
{\rho^p} {\pi} {{{R}}^{{2}}},
```

```math
{I^p}
=
\frac{{1}}{{2}} {\rho^p} {\pi} {{{R}}^{{4}}},
```

with ${R}$ being the radius of the particle.

The velocity vector ${{\Upsilon}_i}$ is used to integrate the particle position in time:

```math
{{{{\Xi}_i}^{n + 1}} = {{{\Xi}_i}^{n}} + \frac{{\Delta {t}}}{{2}} \left( {{{\Upsilon}_i}^{n + 1}} + {{{\Upsilon}_i}^{n}} \right)}.
```

Since ${{\phi}^{n + 1}}$ is a function of ${{{\Xi}_i}^{n + 1}}$, we iterate the equations until the position vector is converged.
Note that ${{{u}_i}^{**}}$ can be used instead of ${{{u}_i}^{n + 1}}$ to evaluate the internal fluid effect, which however makes the scheme less stable.

