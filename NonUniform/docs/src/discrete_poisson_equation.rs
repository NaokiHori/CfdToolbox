use komainu::{MathExpr, doc};

fn one() -> MathExpr {
    MathExpr::from("1")
}

enum Direction {
    X,
    Y,
}

impl AsRef<Direction> for Direction {
    fn as_ref(&self) -> &Self {
        self
    }
}

impl Direction {
    fn to_number(&self) -> &'static str {
        match self {
            Self::X => "1",
            Self::Y => "2",
        }
    }

    fn to_index(&self) -> &'static str {
        match self {
            Self::X => "i",
            Self::Y => "j",
        }
    }
}

enum Offset {
    MinusHalf,
    MinusOne,
    None,
    PlusHalf,
    PlusOne,
}

impl AsRef<Offset> for Offset {
    fn as_ref(&self) -> &Self {
        self
    }
}

impl std::fmt::Display for Offset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::MinusHalf => r"-\frac{1}{2}",
                Self::MinusOne => "-1",
                Self::None => "",
                Self::PlusHalf => r"+\frac{1}{2}",
                Self::PlusOne => "+1",
            }
        )
    }
}

fn linear_system_coefficient<T: AsRef<Direction>>(symbol: &MathExpr, direction: T) -> MathExpr {
    let direction = direction.as_ref();
    MathExpr::from(
        format!(
            r"\left( {}_{} \right)_{{{}}}",
            symbol,
            direction.to_number(),
            direction.to_index(),
        )
        .as_str(),
    )
}

fn at_1d(variable: &MathExpr, offset: Offset) -> MathExpr {
    MathExpr::from(format!(r"{}_{{{}{}}}", variable, Direction::X.to_index(), offset,).as_str())
}

fn at_2d(variable: &MathExpr, offsets: (Offset, Offset)) -> MathExpr {
    MathExpr::from(
        format!(
            r"{}_{{{}{}, {}{}}}",
            variable,
            Direction::X.to_index(),
            offsets.0,
            Direction::Y.to_index(),
            offsets.1
        )
        .as_str(),
    )
}

fn at_x_eigenspace(variable: &MathExpr, offset: Offset) -> MathExpr {
    MathExpr::from(
        format!(
            r"{}_{{\lambda_k, {}{}}}",
            variable,
            Direction::Y.to_index(),
            offset
        )
        .as_str(),
    )
}

fn equal<T: AsRef<MathExpr>, U: AsRef<MathExpr>>(lhs: T, rhs: U) -> MathExpr {
    let lhs = lhs.as_ref();
    let rhs = rhs.as_ref();
    MathExpr::from(&format!("{lhs} = {rhs}"))
}

fn matrix(symbol: &str, indices: (&str, &str)) -> MathExpr {
    MathExpr::from(format!("{}_{{{}{}}}", symbol, &indices.0, &indices.1))
}

pub fn exec(symbols: &crate::Symbols, destination: &str) {
    let coordinate = |direction: &Direction| -> MathExpr {
        MathExpr::from(format!("{}^{}", &symbols.coordinate, direction.to_number()).as_str())
    };
    let scale_factor = |direction: Direction, offset: Offset| -> MathExpr {
        let offset = offset.as_ref();
        let coordinate = coordinate(&direction);
        MathExpr::from(
            format!(
                r"\left( {}_{} \right)_{{{}{}}}",
                &symbols.scale_factor,
                coordinate,
                direction.to_index(),
                offset
            )
            .as_str(),
        )
    };
    let f = MathExpr::from("f");
    let g = MathExpr::from("g");
    let block_diagonal_matrix_symbol = MathExpr::from("B");
    let f_x_eigenspace = MathExpr::from("F");
    let g_x_eigenspace = MathExpr::from("G");
    let f_vector = |index: &str| -> MathExpr { MathExpr::from(format!("{}_{index}", &f)) };
    let g_vector = |index: &str| -> MathExpr { MathExpr::from(format!("{}_{index}", &g)) };
    let f_vector_x_eigenspace =
        |index: &str| -> MathExpr { MathExpr::from(format!("{}_{index}", &f_x_eigenspace)) };
    let g_vector_x_eigenspace =
        |index: &str| -> MathExpr { MathExpr::from(format!("{}_{index}", &g_x_eigenspace)) };
    let lower_symbol = MathExpr::from("l");
    let upper_symbol = MathExpr::from("u");
    let x_eigenvalue = MathExpr::from(r"\lambda_x");
    let x_lower_diagonal = linear_system_coefficient(&lower_symbol, Direction::X);
    let x_upper_diagonal = linear_system_coefficient(&upper_symbol, Direction::X);
    let y_lower_diagonal = linear_system_coefficient(&lower_symbol, Direction::Y);
    let y_upper_diagonal = linear_system_coefficient(&upper_symbol, Direction::Y);
    let laplace_operator_symbol = "L";
    let real_symmetric_matrix_symbol = "S";
    let diagonal_matrix_symbol = r"\Lambda";
    let orthogonal_matrix_symbol = "Q";
    let symmetrization_matrix_symbol = "D";
    let inverse_symmetrization_matrix_symbol = "D^{-1}";
    let laplace_operator =
        |indices: (&str, &str)| -> MathExpr { matrix(laplace_operator_symbol, indices) };
    let real_symmetric_matrix =
        |indices: (&str, &str)| -> MathExpr { matrix(real_symmetric_matrix_symbol, indices) };
    let markdown = doc!(
        include_str!("discrete_poisson_equation.komainu.md"),
        block_diagonal_matrix_symbol = &block_diagonal_matrix_symbol,
        diagonal_matrix = matrix(diagonal_matrix_symbol, ("k", "l")),
        discrete_poisson_equation = equal(
            &x_lower_diagonal * at_2d(&f, (Offset::MinusOne, Offset::None))
                - (&x_lower_diagonal + &x_upper_diagonal) * at_2d(&f, (Offset::None, Offset::None))
                + &x_upper_diagonal * at_2d(&f, (Offset::PlusOne, Offset::None))
                + &y_lower_diagonal * at_2d(&f, (Offset::None, Offset::MinusOne))
                - (&y_lower_diagonal + &y_upper_diagonal) * at_2d(&f, (Offset::None, Offset::None))
                + &y_upper_diagonal * at_2d(&f, (Offset::None, Offset::PlusOne)),
            at_2d(&g, (Offset::None, Offset::None)),
        ),
        discrete_poisson_equation_1d = equal(
            &x_lower_diagonal * at_1d(&f, Offset::MinusOne)
                - (&x_lower_diagonal + &x_upper_diagonal) * at_1d(&f, Offset::None)
                + &x_upper_diagonal * at_1d(&f, Offset::PlusOne),
            at_1d(&g, Offset::None),
        ),
        discrete_poisson_equation_1d_matrix_form =
            equal(laplace_operator(("i", "j")) * &f_vector("j"), g_vector("i")),
        discrete_poisson_equation_1d_matrix_form_decomposed = equal(
            matrix(diagonal_matrix_symbol, ("i", "k"))
                * matrix(orthogonal_matrix_symbol, ("l", "k"))
                * matrix(symmetrization_matrix_symbol, ("l", "j"))
                * &f_vector("j"),
            matrix(orthogonal_matrix_symbol, ("k", "i"))
                * matrix(symmetrization_matrix_symbol, ("k", "j"))
                * &g_vector("j"),
        ),
        discrete_poisson_equation_1d_matrix_form_decomposed_eigenspace = equal(
            matrix(diagonal_matrix_symbol, ("i", "j")) * f_vector_x_eigenspace("j"),
            g_vector_x_eigenspace("i"),
        ),
        discrete_poisson_equation_x_eigenspace = equal(
            &x_eigenvalue * at_x_eigenspace(&f_x_eigenspace, Offset::None)
                + &y_lower_diagonal * at_x_eigenspace(&f_x_eigenspace, Offset::MinusOne)
                - (&y_lower_diagonal + &y_upper_diagonal)
                    * at_x_eigenspace(&f_x_eigenspace, Offset::None)
                + &y_upper_diagonal * at_x_eigenspace(&f_x_eigenspace, Offset::PlusOne),
            at_x_eigenspace(&g_x_eigenspace, Offset::None),
        ),
        f = &f,
        f_x_eigenspace = &f_x_eigenspace,
        g = &g,
        laplace_operator = laplace_operator(("i", "j")),
        orthogonal_matrix = matrix(orthogonal_matrix_symbol, ("i", "k")),
        orthogonal_matrix_symbol = &orthogonal_matrix_symbol,
        real_symmetric_matrix_definition = equal(
            real_symmetric_matrix(("i", "j")),
            matrix(orthogonal_matrix_symbol, ("i", "k"))
                * matrix(diagonal_matrix_symbol, ("k", "l"))
                * matrix(orthogonal_matrix_symbol, ("j", "l"))
        ),
        symmetrization_matrix = matrix(symmetrization_matrix_symbol, ("i", "j")),
        symmetrization_matrix_diagonal = equal(
            matrix(symmetrization_matrix_symbol, ("i", "i")),
            MathExpr::sqrt(scale_factor(Direction::X, Offset::None)),
        ),
        symmetrize_laplace_operator = equal(
            real_symmetric_matrix(("i", "j")),
            matrix(symmetrization_matrix_symbol, ("i", "k"))
                * laplace_operator(("k", "l"))
                * matrix(inverse_symmetrization_matrix_symbol, ("l", "j")),
        ),
        x_eigenvalue = x_eigenvalue,
        x_lower_diagonal_definition = equal(
            &x_lower_diagonal,
            MathExpr::fraction(one(), scale_factor(Direction::X, Offset::None))
                * MathExpr::fraction(one(), scale_factor(Direction::X, Offset::MinusHalf))
        ),
        x_upper_diagonal_definition = equal(
            &x_upper_diagonal,
            MathExpr::fraction(one(), scale_factor(Direction::X, Offset::None))
                * MathExpr::fraction(one(), scale_factor(Direction::X, Offset::PlusHalf))
        ),
        y_lower_diagonal_definition = equal(
            &y_lower_diagonal,
            MathExpr::fraction(one(), scale_factor(Direction::Y, Offset::None))
                * MathExpr::fraction(one(), scale_factor(Direction::Y, Offset::MinusHalf))
        ),
        y_upper_diagonal_definition = equal(
            &y_upper_diagonal,
            MathExpr::fraction(one(), scale_factor(Direction::Y, Offset::None))
                * MathExpr::fraction(one(), scale_factor(Direction::Y, Offset::PlusHalf))
        ),
    );
    std::fs::write(destination, &markdown.content)
        .unwrap_or_else(|_| println!("failed to write file"));
}
