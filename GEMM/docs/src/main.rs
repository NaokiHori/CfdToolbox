use komainu::MathExpr;

mod discrete_poisson_equation;
mod governing_equation;

struct Symbols {
    coordinate: MathExpr,
    scale_factor: MathExpr,
    time: MathExpr,
}

fn main() {
    let symbols = Symbols {
        coordinate: MathExpr::from(r"\xi"),
        scale_factor: MathExpr::from("h"),
        time: MathExpr::from("t"),
    };
    governing_equation::exec(&symbols, "governing_equation/README.md");
    discrete_poisson_equation::exec(&symbols, "discrete_poisson_equation/README.md");
}
