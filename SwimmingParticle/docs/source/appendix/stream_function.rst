
.. _appendix_stream_function:

###############
Stream function
###############

We consider a two-dimensional domain filled with a liquid.
Assuming that the liquid is incompressible and satisfies the Stokes equation:

.. math::

    0_i
    =
    -
    \pder{}{p}{x_i}
    +
    \frac{1}{Re}
    \pder{}{}{x_j}
    \pder{}{u_i}{x_j},

we obtain the biharmonic equation:

.. math::

    \pder{}{}{x_j}
    \pder{}{}{x_j}
    \pder{}{}{x_i}
    \pder{}{\psi}{x_i}
    =
    0,

where :math:`\psi` is the stream function.
The main focus of this section is to compute :math:`\psi` from a set of velocity boundary conditions.

For brevity, we assume the domain is periodic in one of the directions, allowing us to use the Fourier series to expand the stream function.
The other direction is wall-bounded, on which impermeable conditions :math:`u_i n_i = 0` with :math:`n_i` being the wall-normal vector are enforced, as well as the no-slip conditions :math:`u_i t_i = f \left( s \right)` with :math:`t_i` being the wall-tangential vector and :math:`s` being the tangential coordinate.

.. toctree::
    :maxdepth: 1

    ./stream_function/cartesian.rst
    ./stream_function/polar.rst

