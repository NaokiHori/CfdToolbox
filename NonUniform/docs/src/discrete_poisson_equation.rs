use komainu::{Document, Expr};

struct Symbols {
    solution: &'static str,
    right_hand_side: &'static str,
    spectral_solution: &'static str,
    spectral_right_hand_side: &'static str,
}

pub fn exec(
    global_symbols: &crate::Symbols,
    input: &str,
    output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let symbols = Symbols {
        solution: "f",
        right_hand_side: "g",
        spectral_solution: r"\phi",
        spectral_right_hand_side: r"\gamma",
    };
    let template = std::fs::read_to_string(input)?;
    let mut doc = Document::new(template);
    doc.symbol("i", "i");
    doc.symbol("x", "x");
    doc.symbol("dof", "N_x");
    doc.symbol(
        "scale_factor",
        &format!(
            r"{{{}}}_{{{}}}",
            global_symbols.scale_factor, global_symbols.coordinate
        ),
    );
    doc.symbol("lower_diagonal", "l");
    doc.symbol("upper_diagonal", "u");
    doc.symbol("solution", symbols.solution);
    doc.symbol("right_hand_side", symbols.right_hand_side);
    doc.symbol("linear_operator", r"\boldsymbol{L}");
    doc.symbol("real_symmetric_operator", r"\boldsymbol{S}");
    doc.symbol("orthogonal_operator", r"\boldsymbol{Q}");
    doc.symbol("diagonal_operator", r"\boldsymbol{\Lambda}");
    doc.symbol("symmetrization_operator", r"\boldsymbol{D}");
    doc.symbol("solution_vector", &format!(r"\vec{{{}}}", symbols.solution));
    doc.symbol(
        "right_hand_side_vector",
        &format!(r"\vec{{{}}}", symbols.right_hand_side),
    );
    doc.symbol(
        "spectral_solution_vector",
        &format!(r"\vec{{{}}}", symbols.spectral_solution),
    );
    doc.symbol(
        "spectral_right_hand_side_vector",
        &format!(r"\vec{{{}}}", symbols.spectral_right_hand_side),
    );
    doc.symbol("eigenvalue_vector", r"\vec{\lambda}");
    doc.operator("sqrt", 1, |args| Expr::sqrt(&args[0]));
    doc.operator("fraction", 2, |args| Expr::fraction(&args[0], &args[1]));
    doc.operator("index", 2, |args| Expr::subscript(&args[0], &args[1]));
    doc.operator("index_2d", 3, |args| {
        Expr::subscript(&args[0], Expr::list([&args[1], &args[2]]))
    });
    let transpose_symbol = doc.symbol("transpose_symbol", "T");
    doc.operator("transpose", 1, move |args| {
        Expr::superscript(&args[0], &transpose_symbol)
    });
    let inverse_symbol = doc.symbol("inverse_symbol", "-1");
    doc.operator("inverse", 1, move |args| {
        Expr::superscript(&args[0], &inverse_symbol)
    });
    let diagonal_symbol = doc.symbol("diagonal_symbol", r"\text{diag}");
    doc.operator("diagonal", 1, move |args| {
        Expr::function(&diagonal_symbol, args)
    });
    let markdown = doc.render()?;
    std::fs::write(output, markdown)?;
    Ok(())
}
