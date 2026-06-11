# Axisymmetric VOF

![cover image](https://github.com/NaokiHori/CfdToolbox/blob/main/AxisymmetricVOF/cover.jpg)

A Navier-Stokes solver for interfacial two-phase flows under axisymmetric conditions.

## Dependencies

- [C Compiler](https://gcc.gnu.org)
- [GNU Make](https://www.gnu.org/software/make/)
- [OpenMP](https://www.openmp.org) (optional)

## Quick Start

Check out the readme of [the other library](https://github.com/NaokiHori/VerySimpleNSSolver).

## Methodology

We consider two incompressible and immiscible Newtonian liquids (e.g., water and air) which are separated by free and deformable surface.
The liquids are confined inside polar coordinates under axisymmetric condition, or mathematically:

```math
u_2 \equiv 0,
```

```math
\frac{\partial q}{\partial \xi^2} \equiv 0.
```

Hereafter subscripts / superscripts 1, 2, 3 denote the radial, azimuthal, axial directions, respectively.

Their dynamics are governed by the incompressibility constraint:

```math
\frac{1}{J}
\frac{\partial}{\partial \xi^1}
\left(
    \frac{J}{h_{\xi^1}} u_1
\right)
+
\frac{1}{J}
\frac{\partial}{\partial \xi^3}
\left(
    \frac{J}{h_{\xi^3}} u_3
\right)
=
0,
```

the mass balance:

```math
\frac{\partial \rho}{\partial t}
+
\frac{1}{J}
\frac{\partial}{\partial \xi^1}
\left(
  \frac{J}{h_{\xi^1}}
  \rho
  u_1
\right)
+
\frac{1}{J}
\frac{\partial}{\partial \xi^3}
\left(
  \frac{J}{h_{\xi^3}}
  \rho
  u_3
\right)
=
0,
```

and the momentum balance:

```math
\frac{\partial \rho u_1}{\partial t}
=
-
\frac{1}{J}
\frac{\partial}{\partial \xi^1}
\left[
    \left(
        \frac{J}{h_{\xi^1}}
        \rho u_1
    \right)
    u_1
\right]
-
\frac{1}{J}
\frac{\partial}{\partial \xi^3}
\left[
    \left(
        \frac{J}{h_{\xi^3}}
        \rho u_3
    \right)
    u_1
\right]
-
\frac{1}{h_{\xi^1}}
\frac{\partial p}{\partial \xi^1}
+
\frac{1}{J}
\frac{\partial}{\partial \xi^1}
\left(
    \frac{
        J
    }{
        h_{\xi^1}
    }
    \tau_{1 1}
\right)
+
\frac{1}{J}
\frac{\partial}{\partial \xi^3}
\left(
    \frac{
        J
    }{
        h_{\xi^3}
    }
    \tau_{3 1}
\right)
-
\frac{1}{J}
\frac{\partial}{\partial \xi^1}
\left(
    \frac{J}{h_{\xi^1}}
\right)
\tau_{2 2},
```

```math
\frac{\partial \rho u_3}{\partial t}
=
-
\frac{1}{J}
\frac{\partial}{\partial \xi^1}
\left[
    \left(
        \frac{J}{h_{\xi^1}}
        \rho u_1
    \right)
    u_3
\right]
-
\frac{1}{J}
\frac{\partial}{\partial \xi^3}
\left[
    \left(
        \frac{J}{h_{\xi^3}}
        \rho u_3
    \right)
    u_3
\right]
-
\frac{1}{h_{\xi^3}}
\frac{\partial p}{\partial \xi^3}
+
\frac{1}{J}
\frac{\partial}{\partial \xi^1}
\left(
    \frac{
        J
    }{
        h_{\xi^1}
    }
    \tau_{1 3}
\right)
+
\frac{1}{J}
\frac{\partial}{\partial \xi^3}
\left(
    \frac{
        J
    }{
        h_{\xi^3}
    }
    \tau_{3 3}
\right).
```

The whole equations are solved on the general orthogonal coordinate system.

The following numerical techniques are adopted to solve the equations.

- Second-order accurate central-finite-difference scheme for spatial discretization
- Predictor-corrector method (`SMAC` / fractional-step methods) for enforcing incompressibility constraint
- Volume-of-fluid (`THINC`) method for interface capturing 
- Continuum surface-tension force model for describing surface-tension force
- Approximating scheme of the variable-coefficient Poisson equation proposed by Dodd and Ferrante to facilitate the fast Poisson solver
- Energy-consistent discretization for stable integration

For more details, check-out the following pages.

- [Simple TC Solver](https://github.com/NaokiHori/SimpleTCSolver)
- [Simple Bubbly Flow Solver](https://github.com/NaokiHori/SimpleBubblyFlowSolver)

