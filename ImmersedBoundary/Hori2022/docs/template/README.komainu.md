We consider a solid circular object, whose position, translational velocity, and angular velocity of the gravity center are denoted by ${{ vector(particle_center_position, i) }}$, ${{ vector(particle_center_velocity, i) }}$, and ${{ vector(particle_center_angular_velocity, i) }}$, respectively.
Due to the rigid body motion, the velocity vector at any point of the object ${{ vector(particle_surface_position, i) }}$ is described as:

```math
{{
  vector(particle_surface_velocity, i)
  =
  vector(particle_center_velocity, i)
  +
  levi_civita
  *
  vector(particle_center_angular_velocity, j)
  *
  (
    vector(particle_surface_position, k)
    -
    vector(particle_center_position, k)
  )
}}.
```

The fluid velocity field is updated using the Euler-forward scheme, without incorporating the effect of the immersed object:

```math
{{ 
  step_at(vector(fluid_velocity, i), predicted)
  =
  step_at(vector(fluid_velocity, i), current)
  +
  time_step_size
  *
  (
    - 
    step_at(vector(fluid_velocity, j), current)
    *
    partial(
      step_at(vector(fluid_velocity, i), current),
      vector(position, j)
    )
    +
    1 
    / 
    fluid_density
    *
    partial(
      step_at(vector(cauchy_stress_tensor, ij), current),
      vector(position, j)
    )
    +
    vector(external_acceleration, i)
  )
}}.
```

To take into account the effect of the particle by enforcing the desired fluid velocity on the particle surface, fluid and particle momenta are exchanged through:

```math
{{
  vector(immersed_object_acceleration, i)
  =
  heaviside( 1 - step_at(particle_volume_fraction, current) )
  *
  step_at(particle_volume_fraction, current)
  *
  fraction(
    step_at(particle_surface_velocity, current)
    -
    step_at(fluid_velocity, predicted),
    time_step_size
  )
}},
```

where ${{ particle_volume_fraction }}$ is the local volume fraction of the particle.

The fluid velocity is further corrected to account for this effect:

```math
{{
  step_at(vector(fluid_velocity, i), second_predicted)
  =
  step_at(vector(fluid_velocity, i), predicted)
  +
  vector(immersed_object_acceleration, i)
  *
  time_step_size
}}.
```

The same variable ${{ vector(immersed_object_acceleration, i) }}$ is used to evaluate the particle acceleration:

```math
{{
  step_at(vector(particle_center_velocity, i), next)
  =
  step_at(vector(particle_center_velocity, i), current)
  -
  1
  /
  particle_mass
  *
  (
    time_step_size * particle_volume_integral(
      particle_volume,
      fluid_density * vector(immersed_object_acceleration, i),
      particle_volume
    )
    +
    particle_volume_integral(
      particle_volume,
      fluid_density 
        * step_at(particle_volume_fraction, next) 
        * step_at(vector(fluid_velocity, i), next),
      particle_volume
    )
    -
    particle_volume_integral(
      particle_volume,
      fluid_density 
        * step_at(particle_volume_fraction, current) 
        * step_at(vector(fluid_velocity, i), current),
      particle_volume
    )
  )
}},
```

```math
{{
  step_at(vector(particle_center_angular_velocity, i), next)
  =
  step_at(vector(particle_center_angular_velocity, i), current)
  -
  1
  /
  particle_moment_of_inertia
  *
  (
    time_step_size
    *
    particle_volume_integral(
      particle_volume,
      fluid_density
      *
      levi_civita
      *
      (
        vector(particle_surface_position, j)
        -
        vector(particle_center_position, j)
      )
      *
      vector(immersed_object_acceleration, k),
      particle_volume
    )
    +
    particle_volume_integral(
      particle_volume,
      fluid_density
      * 
      step_at(particle_volume_fraction, next)
      *
      levi_civita
      *
      ( 
        vector(particle_surface_position, j)
        -
        vector(particle_center_position, j)
      )
      *
      step_at(vector(fluid_velocity, k), next),
      particle_volume
    )
    -
    particle_volume_integral(
      particle_volume,
      fluid_density
      * 
      step_at(particle_volume_fraction, current)
      *
      levi_civita
      *
      ( 
        vector(particle_surface_position, j)
        -
        vector(particle_center_position, j)
      )
      *
      step_at(vector(fluid_velocity, k), current),
      particle_volume
    )
  )
}},
```

where ${{ particle_mass }}$ and ${{ particle_moment_of_inertia }}$ are the mass and the moment of inertia of the particle:

```math
{{ particle_mass = particle_density * pi * particle_radius^2 }},
```

```math
{{ particle_moment_of_inertia = 1 / 2 * particle_density * pi * particle_radius^4 }},
```

with ${{ particle_radius }}$ being the radius of the particle.

The velocity vector ${{ vector(particle_center_velocity, i) }}$ is used to integrate the particle position in time:

```math
{{
  step_at(vector(particle_center_position, i), next)
  =
  step_at(vector(particle_center_position, i), current)
  +
  time_step_size
  /
  2
  *
  (
    step_at(vector(particle_center_velocity, i), next)
    +
    step_at(vector(particle_center_velocity, i), current)
  )
}}.
```

Since ${{ step_at(particle_volume_fraction, next) }}$ is a function of ${{ step_at(vector(particle_center_position, i), next) }}$, we iterate the equations until the position vector is converged.
Note that ${{ step_at(vector(fluid_velocity, i), second_predicted) }}$ can be used instead of ${{ step_at(vector(fluid_velocity, i), next) }}$ to evaluate the internal fluid effect, which however makes the scheme less stable.

