use komainu::{MathExpr, doc};

enum Step {
    Current,
    Next,
    Predicted,
}

impl std::fmt::Display for Step {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        let step = match self {
            Self::Current => "n",
            Self::Next => "n + 1",
            Self::Predicted => "*",
        };
        let step = MathExpr::from(step);
        write!(f, "{step}")
    }
}

fn zero() -> MathExpr {
    MathExpr::from("0")
}

fn one() -> MathExpr {
    MathExpr::from("1")
}

fn at<T: AsRef<MathExpr>>(target: T, step: Step) -> MathExpr {
    let target = target.as_ref();
    MathExpr::from(format!("{target}^{step}").as_str())
}

fn velocity(index: &str) -> MathExpr {
    let symbol = MathExpr::from("u");
    let index = MathExpr::from(index);
    MathExpr::from(format!("{symbol}_{index}").as_str())
}

fn equal<T: AsRef<MathExpr>, U: AsRef<MathExpr>>(lhs: T, rhs: U) -> MathExpr {
    let lhs = lhs.as_ref();
    let rhs = rhs.as_ref();
    MathExpr::from(&format!("{lhs} = {rhs}"))
}

pub fn exec(symbols: &crate::Symbols, destination: &str) {
    let coordinate = |index: &str| -> MathExpr {
        let index = MathExpr::from(index);
        MathExpr::from(format!("{}^{}", &symbols.coordinate, index).as_str())
    };
    let scale_factor = |index: &str| -> MathExpr {
        let coordinate = coordinate(index);
        MathExpr::from(format!("{}_{}", &symbols.scale_factor, coordinate).as_str())
    };
    let pressure = MathExpr::from("p");
    let scalar_potential = MathExpr::from(r"\psi");
    let jacobian_determinant = MathExpr::from("J");
    let cartesian = MathExpr::from("x_i");
    let time_step_size = MathExpr::from(format!(r"\Delta {}", &symbols.time));
    let coordinate_1 = coordinate("1");
    let coordinate_2 = coordinate("2");
    let scale_factor_1 = scale_factor("1");
    let scale_factor_2 = scale_factor("2");
    //
    let markdown = doc!(
        include_str!("governing_equation.komainu.md"),
        coordinate_1 = coordinate_1,
        coordinate_2 = coordinate_2,
        correct_velocity = equal(
            at(velocity("i"), Step::Next),
            at(velocity("i"), Step::Predicted)
                - &time_step_size
                    * MathExpr::fraction(one(), scale_factor("i"))
                    * MathExpr::partial_derivative(&scalar_potential, coordinate("i"))
        ),
        incompressibility_constraint = equal(
            MathExpr::fraction(one(), &jacobian_determinant)
                * MathExpr::partial_derivative_operator(
                    MathExpr::fraction(&jacobian_determinant, scale_factor("j")) * velocity("j"),
                    coordinate("j")
                ),
            zero()
        ),
        jacobian_determinant = &jacobian_determinant,
        jacobian_determinant_definition = format!(
            r"{} \equiv \Pi_i {}",
            jacobian_determinant,
            scale_factor("i")
        ),
        momentum_balance = equal(
            MathExpr::partial_derivative(velocity("i"), &symbols.time),
            -MathExpr::fraction(velocity("j"), scale_factor("j"))
                * MathExpr::partial_derivative(velocity("i"), coordinate("j"))
                - MathExpr::fraction(one(), scale_factor("i"))
                    * MathExpr::partial_derivative(&pressure, coordinate("i"))
                + MathExpr::fraction(one(), &jacobian_determinant)
                    * MathExpr::partial_derivative_operator(
                        MathExpr::fraction(&jacobian_determinant, scale_factor("j"))
                            * MathExpr::fraction(one(), scale_factor("j"))
                            * MathExpr::partial_derivative(velocity("i"), coordinate("j")),
                        coordinate("i")
                    )
        ),
        poisson_equation = equal(
            MathExpr::fraction(one(), &jacobian_determinant)
                * MathExpr::partial_derivative_operator(
                    MathExpr::fraction(&jacobian_determinant, scale_factor("j"))
                        * MathExpr::fraction(one(), scale_factor("j"))
                        * MathExpr::partial_derivative(&scalar_potential, coordinate("j")),
                    coordinate("j")
                ),
            MathExpr::fraction(one(), &time_step_size)
                * MathExpr::fraction(one(), &jacobian_determinant)
                * MathExpr::partial_derivative_operator(
                    MathExpr::fraction(&jacobian_determinant, scale_factor("j"))
                        * at(velocity("j"), Step::Predicted),
                    coordinate("j")
                ),
        ),
        predict_velocity = equal(
            at(velocity("i"), Step::Predicted),
            at(velocity("i"), Step::Current)
                + MathExpr::integral(
                    Some(at(&symbols.time, Step::Current)),
                    Some(at(&symbols.time, Step::Next)),
                    MathExpr::from(r"\left( RHS \right)_i"),
                    &symbols.time,
                )
        ),
        scalar_potential = &scalar_potential,
        scale_factor = scale_factor("i"),
        scale_factor_1 = &scale_factor_1,
        scale_factor_2 = scale_factor("2"),
        scale_factor_definition = format!(
            r"{} \equiv {}",
            scale_factor("i"),
            MathExpr::partial_derivative(&cartesian, coordinate("i")),
        ),
        simplified_poisson_equation = equal(
            MathExpr::fraction(one(), &scale_factor_1)
                * MathExpr::partial_derivative_operator(
                    MathExpr::fraction(one(), &scale_factor_1)
                        * MathExpr::partial_derivative(&scalar_potential, &coordinate_1),
                    &coordinate_1
                )
                + MathExpr::fraction(one(), &scale_factor_2)
                    * MathExpr::partial_derivative_operator(
                        MathExpr::fraction(one(), &scale_factor_2)
                            * MathExpr::partial_derivative(&scalar_potential, &coordinate_2),
                        &coordinate_2
                    ),
            MathExpr::fraction(one(), &time_step_size)
                * (MathExpr::fraction(one(), &scale_factor_1)
                    * MathExpr::partial_derivative(
                        at(velocity("1"), Step::Predicted),
                        coordinate_1
                    )
                    + MathExpr::fraction(one(), &scale_factor_2)
                        * MathExpr::partial_derivative(
                            at(velocity("2"), Step::Predicted),
                            coordinate_2
                        )),
        ),
    );
    std::fs::write(destination, &markdown.content)
        .unwrap_or_else(|_| println!("failed to write file"));
}
