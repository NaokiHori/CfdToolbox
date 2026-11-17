// mod discrete_poisson_equation;
mod governing_equation;

struct Symbols {
    coordinate: String,
    scale_factor: String,
    time: String,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let symbols = Symbols {
        coordinate: String::from(r"\xi"),
        scale_factor: String::from("h"),
        time: String::from("t"),
    };
    governing_equation::exec(&symbols, "governing_equation.md")?;
    // discrete_poisson_equation::exec(&symbols, "discrete_poisson_equation/README.md");
    Ok(())
}
