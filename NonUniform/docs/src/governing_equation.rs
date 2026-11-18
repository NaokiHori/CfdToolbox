use komainu::{Document, Expr};

fn coordinate(symbol: &Expr, index: &Expr) -> Expr {
    Expr::superscript(symbol, index)
}

pub fn exec(
    symbols: &crate::Symbols,
    input: &str,
    output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let template = std::fs::read_to_string(input)?;
    let mut doc = Document::new(template);
    doc.symbol("i", r"i");
    doc.symbol("j", r"j");
    doc.symbol("time", &symbols.time);
    doc.symbol("jacobian_determinant", r"J");
    doc.symbol("reynolds", r"Re");
    doc.symbol("pressure", r"p");
    doc.symbol("scalar_potential", r"\psi");
    doc.symbol("right_hand_side", r"(RHS)");
    doc.symbol("predicted", r"*");
    doc.symbol("current", r"n");
    doc.symbol("next", r"n + 1");
    let cartesian_symbol = doc.symbol("cartesian_symbol", r"x");
    let coordinate_symbol = doc.symbol("coordinate_symbol", &symbols.coordinate);
    let velocity_symbol = doc.symbol("velocity_symbol", r"u");
    let scale_factor_symbol = doc.symbol("scale_factor_symbol", &symbols.scale_factor);
    doc.operator("vector", 2, |args| Expr::subscript(&args[0], &args[1]));
    doc.operator("cartesian", 1, {
        let cartesian_symbol = cartesian_symbol.clone();
        move |args| Expr::subscript(&cartesian_symbol, &args[0])
    });
    doc.operator("coordinate", 1, {
        let coordinate_symbol = coordinate_symbol.clone();
        move |args| coordinate(&coordinate_symbol, &args[0])
    });
    doc.operator("velocity", 1, {
        let velocity_symbol = velocity_symbol.clone();
        move |args| Expr::subscript(&velocity_symbol, &args[0])
    });
    doc.operator("scale_factor", 1, {
        let coordinate_symbol = coordinate_symbol.clone();
        move |args| {
            Expr::subscript(
                &scale_factor_symbol,
                coordinate(&coordinate_symbol, &args[0]),
            )
        }
    });
    doc.operator("partial", 2, |args| {
        Expr::partial_derivative(&args[0], &args[1])
    });
    doc.operator("partial_operator", 2, |args| {
        Expr::partial_derivative_operator(&args[0], &args[1])
    });
    doc.operator("fraction", 2, |args| Expr::fraction(&args[0], &args[1]));
    doc.operator("product_over", 2, |args| {
        Expr::product_over(&args[0], &args[1])
    });
    doc.operator("step_at", 2, |args| Expr::superscript(&args[0], &args[1]));
    doc.operator("definite_integral", 4, |args| {
        Expr::definite_integral(&args[0], &args[1], &args[2], &args[3])
    });
    let markdown = doc.render()?;
    std::fs::write(output, markdown)?;
    Ok(())
}
