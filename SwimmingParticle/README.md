# Swimming particle

![cover image](https://github.com/NaokiHori/CfdToolbox/blob/main/SwimmingParticle/cover.gif)

## Overview

A finite-difference-based numerical simulator for diffusophoretic, self-propelling particles.

## Quick Start

Build and run:

```console
make output
make all
./a.out
```

Check the result:

```console
python visualize.py
```

## Parameters

There are four parameters which can be changed:

- Péclet number
- Radius of the outer circle
- Number of radial and azimuthal grid points 

They are all defined in `src/main.c`.
**Note that the number of azimuthal grids are to be a power of 2.**

Re-compile the source (`make all`) when modified.

## Documentation

Governing equation and its numerical treatment are briefly documented in `docs`.

For now the documentation is not hosted, so build it by yourself if needed.
Follow [README.md](./docs/README.md).

