# GEMM

![cover.jpg](./cover.jpg)

CFD solver incorporating the general matrix-matrix multiplication technique, enabling the use of stretched grids in multiple directions.

## Dependency

For now eigendecomposition is out-sourced to NumPy:

- Python3 (NumPy)

## Quick Start

Set up a grid:

```console
python3 setup-grid.py
```

The resulting grid is saved to `grid/`.
Build solver and run:

```console
make output
make all
./a.out
```

## Caveat

The objective is to understand how things go.
Performance is out-of consideration.

## Reference

- [Costa et al., arXiv, 2026](https://arxiv.org/abs/2603.09528)

