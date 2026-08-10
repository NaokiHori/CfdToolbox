# Governing equation

The equations are expressed in the general orthogonal coordinate system:

```math
{\frac{{1}}{{J}} \frac{\partial}{\partial {{\xi}^{j}}} \left( \frac{{J}}{{{h}_{{\xi}^{j}}}} {{u}_{j}} \right) = {0}},
```

```math
{\frac{\partial {{u}_{i}}}{\partial {t}} = - \frac{{{u}_{j}}}{{{h}_{{\xi}^{j}}}} \frac{\partial {{u}_{i}}}{\partial {{\xi}^{j}}} - \frac{{1}}{{{h}_{{\xi}^{i}}}} \frac{\partial {p}}{\partial {{\xi}^{i}}} + \frac{{1}}{{J}} \frac{\partial}{\partial {{\xi}^{i}}} \left( \frac{{J}}{{{h}_{{\xi}^{j}}}} \frac{{1}}{{{h}_{{\xi}^{j}}}} \frac{\partial {{u}_{i}}}{\partial {{\xi}^{j}}} \right)},
```

where ${{h}_{{\xi}^{i}}}$ is the scale factor:

```math
{{h}_{{\xi}^{i}}} \equiv \frac{\partial {x_i}}{\partial {{\xi}^{i}}}
```

and ${J}$ is the Jacobian determinant:

```math
{J} \equiv \Pi_i {{h}_{{\xi}^{i}}}.
```

Note that summation is taken with respect to the subscript $j$.
See [the other documentation](https://naokihori.github.io/SimpleTCSolver/index.html) for details.

The equations are integrated by means of the predictor-corrector approach: the velocity field is predicted:

```math
{{{{u}_{i}}^{*}} = {{{u}_{i}}^{n}} + {\int_{{{t}^{n}}}^{{{t}^{n + 1}}} {\left( RHS \right)_i} d{t}}},
```

which is corrected to enforce the incompressibility:

```math
{{{{u}_{i}}^{n + 1}} = {{{u}_{i}}^{*}} - {\Delta {t}} \frac{{1}}{{{h}_{{\xi}^{i}}}} \frac{\partial {\psi}}{\partial {{\xi}^{i}}}},
```

where ${\psi}$ is scalar field obtained by solving the Poisson equation:

```math
{\frac{{1}}{{J}} \frac{\partial}{\partial {{\xi}^{j}}} \left( \frac{{J}}{{{h}_{{\xi}^{j}}}} \frac{{1}}{{{h}_{{\xi}^{j}}}} \frac{\partial {\psi}}{\partial {{\xi}^{j}}} \right) = \frac{{1}}{{\Delta {t}}} \frac{{1}}{{J}} \frac{\partial}{\partial {{\xi}^{j}}} \left( \frac{{J}}{{{h}_{{\xi}^{j}}}} {{{u}_{j}}^{*}} \right)}.
```

Now, we write down this equation explicitly in a two-dimensional Cartesian coordinates.
The metrics have the following relations:

```math
{{h}_{{\xi}^{1}}}
=
{{h}_{{\xi}^{1}}} \left( {{\xi}^{1}} \right),
```

```math
{{h}_{{\xi}^{2}}}
=
{{h}_{{\xi}^{2}}} \left( {{\xi}^{2}} \right),
```

```math
{J}
=
{{h}_{{\xi}^{1}}} {{h}_{{\xi}^{2}}}.
```

Since the scale factor in the first direction is independent in the second direction (and vice versa), the Poisson equation can be simplified as:

```math
{\frac{{1}}{{{h}_{{\xi}^{1}}}} \frac{\partial}{\partial {{\xi}^{1}}} \left( \frac{{1}}{{{h}_{{\xi}^{1}}}} \frac{\partial {\psi}}{\partial {{\xi}^{1}}} \right) + \frac{{1}}{{{h}_{{\xi}^{2}}}} \frac{\partial}{\partial {{\xi}^{2}}} \left( \frac{{1}}{{{h}_{{\xi}^{2}}}} \frac{\partial {\psi}}{\partial {{\xi}^{2}}} \right) = \frac{{1}}{{\Delta {t}}} \left( \frac{{1}}{{{h}_{{\xi}^{1}}}} \frac{\partial {{{u}_{1}}^{*}}}{\partial {{\xi}^{1}}} + \frac{{1}}{{{h}_{{\xi}^{2}}}} \frac{\partial {{{u}_{2}}^{*}}}{\partial {{\xi}^{2}}} \right)}.
```

