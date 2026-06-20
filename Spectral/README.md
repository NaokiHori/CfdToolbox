# Spectral

[![CoverImage](https://github.com/NaokiHori/CfdToolbox/blob/main/Spectral/cover.jpg)](https://youtu.be/U-Wwif6s9yM)

## Overview

Incompressible Navier-Stokes solver by means of spectral method.

## Equation

Using sinusoidal functions as basis functions, we consider the Navier-Stokes equation in the spectral domain:

```math
\frac{d {U_{ij}}}{d {t}} = P_{ij} - \frac{{\xi}_{i}}{{{\xi}_{i}}^{2} + {{\eta}_{j}}^{2}} ({\xi}_{i} P_{ij} + {\eta}_{j} Q_{ij}) - \frac{1}{Re} ({{\xi}_{i}}^{2} + {{\eta}_{j}}^{2}) U_{ij},
```

```math
\frac{d {V_{ij}}}{d {t}} = Q_{ij} - \frac{{\eta}_{j}}{{{\xi}_{i}}^{2} + {{\eta}_{j}}^{2}} ({\xi}_{i} P_{ij} + {\eta}_{j} Q_{ij}) - \frac{1}{Re} ({{\xi}_{i}}^{2} + {{\eta}_{j}}^{2}) V_{ij},
```

where `i = 0, 1, ..., Nx / 2 + 1` and `j = 0, 1, ..., Ny - 1` are integers used to distinguish different wavenumbers.
The wavenumbers in x direction are defined as:

```math
\xi_i = 2 \pi i / L_x,
```

while in y:

```math
\eta_j = 2 \pi j / L_y
```

for `j < Ny / 2 + 1` and 

```math
\eta_j = 2 \pi (j - N_y) / L_y
```

for `Ny / 2 < j < Ny - 1`.

`P` and `Q` are non-linear terms which are defined as:

```math
P_{ij} = - I {\xi}_{i} \sum_{i_0 + i_1 = i} \sum_{j_0 + j_1 = j} U_{i_0 j_0} U_{i_1 j_1} - I {\eta}_{j} \sum_{i_0 + i_1 = i} \sum_{j_0 + j_1 = j} V_{i_0 j_0} U_{i_1 j_1},
```

```math
Q_{ij} = - I {\xi}_{i} \sum_{i_0 + i_1 = i} \sum_{j_0 + j_1 = j} U_{i_0 j_0} V_{i_1 j_1} - I {\eta}_{j} \sum_{i_0 + i_1 = i} \sum_{j_0 + j_1 = j} V_{i_0 j_0} V_{i_1 j_1}.
```

Note that the convolutions are evaluated in the physical space, with 2/3 dealiasing.

## Time evolution

The equations for each wavenumber have the form:

```math
\frac{d p}{d t} + C p = q,
```

which can be solved by the integrating-factor technique:

```math
\frac{d}{d t} [p e^{C (t - t_n)}] = q e^{C (t - t_n)}.
```

This is used along with the low-storage Runge-Kutta method proposed by Williamson 1980; namely, for each Runge-Kutta step, we update two buffers `p` and `g`:

```math
g^* \leftarrow \alpha^k g + q^k \Delta t,
```

```math
p^* \leftarrow p^k + \beta^k g^*,
```

```math
p^{k + 1} \leftarrow p^* e^{- C \gamma^k \Delta t},
g^{k + 1} \leftarrow g^* e^{- C \gamma^k \Delta t},
```

where

```math
p^0 \equiv p (t_n),
p^3 \equiv p (t_{n + 1}).
```

The coefficients are adjusted to achieve the third-order accuracy:

```math
\alpha^0 = 0,
\alpha^1 = - 5 / 9,
\alpha^2 = - 153 / 128,
```

```math
\beta^0 = 1 / 3,
\beta^1 = 15 / 16,
\beta^2 = 8 / 15,
```

```math
\gamma^0 = 1 / 3,
\gamma^1 = 5 / 12,
\gamma^2 = 1 / 4.
```

## See also

See [the other repository](https://github.com/NaokiHori/SpectralNSSolver1) for detailed derivation of the governing equations.

## Reference

- Canuto et al., *Spectral Methods - Fundamentals in Single Domains*, Springer
- Canuto et al., *Spectral Methods - Evolution to Complex Geometries and Applications to Fluid Dynamics*, Springer
- Williamson, J. Comput. Phys. (35), 1980

