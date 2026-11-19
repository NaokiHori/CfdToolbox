mod discrete_poisson_equation;
mod governing_equation;

struct Symbols {
    coordinate: &'static str,
    scale_factor: &'static str,
    time: &'static str,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let symbols = Symbols {
        coordinate: r"\xi",
        scale_factor: "h",
        time: "t",
    };
    governing_equation::exec(
        &symbols,
        "template/governing_equation.komainu.md",
        "output/governing_equation.md",
    )?;
    discrete_poisson_equation::exec(
        &symbols,
        "template/discrete_poisson_equation.komainu.md",
        "output/discrete_poisson_equation.md",
    )?;
    Ok(())
}
