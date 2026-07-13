use komainu::{MathExpr, doc};

fn one() -> MathExpr {
    MathExpr::from("1")
}

// TODO: flexible index
fn vector<T: AsRef<MathExpr>>(symbol: T, index: &str) -> MathExpr {
    let symbol = symbol.as_ref();
    MathExpr::from(format!("{symbol}_{index}").as_str())
}

fn at<T: AsRef<MathExpr>>(variable: T, step: Step) -> MathExpr {
    let variable = variable.as_ref();
    MathExpr::from(format!("{variable}^{step}").as_str())
}

fn volume_integral<T: AsRef<MathExpr>>(integrand: T) -> MathExpr {
    MathExpr::integral(
        None::<MathExpr>,
        None::<MathExpr>,
        integrand,
        MathExpr::from("V"),
    )
}

fn outer_product<T: AsRef<MathExpr>, U: AsRef<MathExpr>>(
    levi_civita: &str,
    vector0: T,
    vector1: U,
) -> MathExpr {
    MathExpr::from(&format!(r"\epsilon_{{{levi_civita}}}")) * vector0.as_ref() * vector1.as_ref()
}

fn heaviside(x: &MathExpr) -> MathExpr {
    MathExpr::from(&format!(r"H \left( {x} \right)"))
}

fn equal<T: AsRef<MathExpr>, U: AsRef<MathExpr>>(lhs: T, rhs: U) -> MathExpr {
    let lhs = lhs.as_ref();
    let rhs = rhs.as_ref();
    MathExpr::from(&format!("{lhs} = {rhs}"))
}

struct Fluid {
    density: MathExpr,
    velocity_vector: MathExpr,
}

struct Particle {
    density: MathExpr,
    center_position_vector: MathExpr,
    center_velocity_vector: MathExpr,
    center_angular_velocity_vector: MathExpr,
    radius: MathExpr,
    surface_position_vector: MathExpr,
    surface_velocity_vector: MathExpr,
    volume_fraction: MathExpr,
    mass: MathExpr,
    moment_of_inertia: MathExpr,
}

enum Step {
    Current,
    FirstPredicted,
    Next,
    SecondPredicted,
}

impl std::fmt::Display for Step {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        let step = match self {
            Self::Current => "n",
            Self::FirstPredicted => "*",
            Self::Next => "n + 1",
            Self::SecondPredicted => "**",
        };
        let step = MathExpr::from(step);
        write!(f, "{step}")
    }
}

fn main() {
    let pi = MathExpr::from(r"\pi");
    let time = MathExpr::from("t");
    let coordinate_symbol = MathExpr::from("x");
    let acceleration_symbol = MathExpr::from("a");
    let velocity_symbol = MathExpr::from("u");
    let gravity_symbol = MathExpr::from("g");
    let cauchy_stress_tensor = MathExpr::from(r"\sigma_{ij}");
    let fluid = Fluid {
        density: MathExpr::from(r"\rho^f"),
        velocity_vector: vector(&velocity_symbol, "i"),
    };
    let particle_center_angular_velocity_vector_symbol = MathExpr::from(r"\Omega");
    let particle_center_position_vector_symbol = MathExpr::from(r"\Xi");
    let particle_surface_position_vector_symbol = MathExpr::from("X");
    let particle = Particle {
        density: MathExpr::from(r"\rho^p"),
        center_position_vector: vector(&particle_center_position_vector_symbol, "i"),
        center_velocity_vector: vector(MathExpr::from(r"\Upsilon"), "i"),
        center_angular_velocity_vector: vector(
            &particle_center_angular_velocity_vector_symbol,
            "i",
        ),
        mass: MathExpr::from("m^p"),
        moment_of_inertia: MathExpr::from("I^p"),
        radius: MathExpr::from("R"),
        surface_position_vector: vector(&particle_surface_position_vector_symbol, "i"),
        surface_velocity_vector: vector(MathExpr::from("U"), "i"),
        volume_fraction: MathExpr::from(r"\phi"),
    };
    let gravity = vector(&gravity_symbol, "i");
    let immersed_object_acceleration = |index: &str| -> MathExpr {
        MathExpr::from(&format!("{}^{{IBM}}", vector(&acceleration_symbol, index)))
    };
    let time_step_size = MathExpr::from(&format!(r"\Delta {time}"));
    let fluid_velocity_vector_current_step =
        |index: &str| -> MathExpr { at(vector(&velocity_symbol, index), Step::Current) };
    let fluid_velocity_vector_first_predicted = at(&fluid.velocity_vector, Step::FirstPredicted);
    let fluid_velocity_vector_second_predicted =
        |index: &str| -> MathExpr { at(vector(&velocity_symbol, index), Step::SecondPredicted) };
    let fluid_velocity_vector_next_step =
        |index: &str| -> MathExpr { at(vector(&velocity_symbol, index), Step::Next) };
    let particle_center_position_vector_next_step =
        at(&particle.center_position_vector, Step::Next);
    let particle_center_velocity_vector_current_step =
        at(&particle.center_velocity_vector, Step::Current);
    let particle_center_velocity_vector_next_step =
        at(&particle.center_velocity_vector, Step::Next);
    let particle_volume_fraction_current_step = at(&particle.volume_fraction, Step::Current);
    let particle_volume_fraction_next_step = at(&particle.volume_fraction, Step::Next);
    //
    let markdown = doc!(
        include_str!("main.komainu.md"),
        fluid_velocity_vector_current_step = fluid_velocity_vector_current_step("i"),
        fluid_velocity_vector_first_predicted = fluid_velocity_vector_first_predicted,
        fluid_velocity_vector_next_step = fluid_velocity_vector_next_step("i"),
        fluid_velocity_vector_second_predicted = fluid_velocity_vector_second_predicted("i"),
        fluid_velocity_vector_update_right_hand_side =
            -at(vector(&velocity_symbol, "j"), Step::Current)
                * MathExpr::partial_derivative(
                    fluid_velocity_vector_current_step("i"),
                    vector(&coordinate_symbol, "j")
                )
                + MathExpr::fraction(one(), &fluid.density)
                    * MathExpr::partial_derivative(
                        at(&cauchy_stress_tensor, Step::Current),
                        vector(&coordinate_symbol, "j")
                    )
                + &gravity,
        immersed_object_acceleration = immersed_object_acceleration("i"),
        immersed_object_acceleration_as_velocity_vector_difference =
            heaviside(&(one() - &particle_volume_fraction_current_step))
                * &particle_volume_fraction_current_step
                * MathExpr::fraction(
                    &(at(&particle.surface_velocity_vector, Step::Current)
                        - fluid_velocity_vector_first_predicted),
                    &time_step_size
                ),
        particle_center_angular_velocity_vector = &particle.center_angular_velocity_vector,
        particle_center_angular_velocity_vector_update = equal(
            at(&particle.center_angular_velocity_vector, Step::Next),
            at(&particle.center_angular_velocity_vector, Step::Current)
                - MathExpr::fraction(one(), &particle.moment_of_inertia)
                    * (&time_step_size
                        * volume_integral(
                            &fluid.density
                                * outer_product(
                                    "ijk",
                                    vector(&particle_surface_position_vector_symbol, "j")
                                        - vector(&particle_center_position_vector_symbol, "j"),
                                    immersed_object_acceleration("k")
                                )
                        )
                        + volume_integral(
                            &fluid.density
                                * &particle_volume_fraction_next_step
                                * outer_product(
                                    "ijk",
                                    vector(&particle_surface_position_vector_symbol, "j")
                                        - vector(&particle_center_position_vector_symbol, "j"),
                                    fluid_velocity_vector_next_step("k")
                                )
                        )
                        - volume_integral(
                            &fluid.density
                                * &particle_volume_fraction_next_step
                                * outer_product(
                                    "ijk",
                                    vector(&particle_surface_position_vector_symbol, "j")
                                        - vector(&particle_center_position_vector_symbol, "j"),
                                    fluid_velocity_vector_current_step("k")
                                )
                        ))
        ),
        particle_center_position_vector = &particle.center_position_vector,
        particle_center_position_vector_next_step = particle_center_position_vector_next_step,
        particle_center_position_vector_update = equal(
            &particle_center_position_vector_next_step,
            &at(&particle.center_position_vector, Step::Current)
                + &(MathExpr::fraction(&time_step_size, MathExpr::from("2"))
                    * (&particle_center_velocity_vector_next_step
                        + &particle_center_velocity_vector_current_step))
        ),
        particle_center_velocity_vector_update = equal(
            &particle_center_velocity_vector_next_step,
            &particle_center_velocity_vector_current_step
                - MathExpr::fraction(one(), &particle.mass)
                    * (&time_step_size
                        * volume_integral(&fluid.density * &immersed_object_acceleration("i"))
                        + volume_integral(
                            &fluid.density
                                * &particle_volume_fraction_next_step
                                * fluid_velocity_vector_next_step("i")
                        )
                        - volume_integral(
                            &fluid.density
                                * &particle_volume_fraction_current_step
                                * fluid_velocity_vector_current_step("i")
                        ))
        ),
        particle_center_velocity_vector = &particle.center_velocity_vector,
        particle_mass = &particle.mass,
        particle_mass_formula = &particle.density * &pi * particle.radius.pow(MathExpr::from("2")),
        particle_moment_of_inertia = &particle.moment_of_inertia,
        particle_moment_of_inertia_formula = MathExpr::fraction(one(), MathExpr::from("2"))
            * &particle.density
            * &pi
            * particle.radius.pow(MathExpr::from("4")),
        particle_radius = &particle.radius,
        particle_surface_position_vector = &particle.surface_position_vector,
        particle_surface_velocity_vector = &particle.surface_velocity_vector,
        particle_surface_velocity_vector_rotational_contribution = outer_product(
            "ijk",
            vector(&particle_center_angular_velocity_vector_symbol, "j"),
            vector(&particle_surface_position_vector_symbol, "k")
                - vector(&particle_center_position_vector_symbol, "k"),
        ),
        particle_volume_fraction = &particle.volume_fraction,
        particle_volume_fraction_next_step = &particle_volume_fraction_next_step,
        time_step_size = time_step_size,
    );
    std::fs::write("README.md", &markdown.content)
        .unwrap_or_else(|_| println!("failed to write file"));
}
