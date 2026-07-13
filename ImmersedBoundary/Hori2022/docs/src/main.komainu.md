We consider a solid circular object, whose position, translational velocity, and angular velocity of the gravity center are denoted by ${particle_center_position_vector}$, ${particle_center_velocity_vector}$, and ${particle_center_angular_velocity_vector}$, respectively.
Due to the rigid body motion, the velocity vector at any point of the object ${particle_surface_position_vector}$ is described as:

```math
{particle_surface_velocity_vector}
=
{particle_center_velocity_vector}
+
{particle_surface_velocity_vector_rotational_contribution}.
```

The fluid velocity field is updated using the Euler-forward scheme, without incorporating the effect of the immersed object:

```math
{fluid_velocity_vector_first_predicted}
=
{fluid_velocity_vector_current_step}
+
{time_step_size}
\left(
    {fluid_velocity_vector_update_right_hand_side}
\right).
```

To take into account the effect of the particle by enforcing the desired fluid velocity on the particle surface, fluid and particle momenta are exchanged through:

```math
{immersed_object_acceleration}
=
{immersed_object_acceleration_as_velocity_vector_difference},
```

where ${particle_volume_fraction}$ is the local volume fraction of the particle.

The fluid velocity is further corrected to account for this effect:

```math
{fluid_velocity_vector_second_predicted}
=
{fluid_velocity_vector_first_predicted}
+
{immersed_object_acceleration}
{time_step_size}.
```

The same variable ${immersed_object_acceleration}$ is used to evaluate the particle acceleration:

```math
{particle_center_velocity_vector_update},
```

```math
{particle_center_angular_velocity_vector_update},
```

where ${particle_mass}$ and ${particle_moment_of_inertia}$ are the mass and the moment of inertia of the particle:

```math
{particle_mass}
=
{particle_mass_formula},
```

```math
{particle_moment_of_inertia}
=
{particle_moment_of_inertia_formula},
```

with ${particle_radius}$ being the radius of the particle.

The velocity vector ${particle_center_velocity_vector}$ is used to integrate the particle position in time:

```math
{particle_center_position_vector_update}.
```

Since ${particle_volume_fraction_next_step}$ is a function of ${particle_center_position_vector_next_step}$, we iterate the equations until the position vector is converged.
Note that ${fluid_velocity_vector_second_predicted}$ can be used instead of ${fluid_velocity_vector_next_step}$ to evaluate the internal fluid effect, which however makes the scheme less stable.

