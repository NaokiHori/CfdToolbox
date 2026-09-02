# Non-Uniform

![cover.jpg](./cover.jpg)

CFD solver incorporating the general matrix-matrix multiplication technique, enabling the use of stretched (non-uniform) grids in multiple directions.

## Quick Start

Build solver and run:

```console
make output
make all
./a.out
```

Domain (lengths, degree of freedoms, grid stretching, among others) is configured in `src/domain.c`.

## Features

- Active scalar coupling (temperature and salinity: double-diffusive convections)
- Three-step Runge-Kutta time integrator
- Implicit diffusive treatment (using factorization)

## Documentation

Visit [docs](./docs).
Important feature is briefly documented.

## Caveat

The objective is to understand how the Poisson equation is handled on non-uniform grid spacings.
Performance is out-of consideration.

## Reference

- [Costa et al., arXiv, 2026](https://arxiv.org/abs/2603.09528)

