# Hori 2022

Immersed boundary method by Hori 2022, which is minified.

## Quick Start

Build and run:

```console
make output
make all
./a.out
```

## Caveats

There are several simplifications:

- Euler forward instead of Adams-Bashforth
- Single particle
- No collision model
- Particle updating strategy (sub iteration only involving internal momentum integration)

## Reference

- Hori et al., _Comput. Fluids_ (**236**), 2022
