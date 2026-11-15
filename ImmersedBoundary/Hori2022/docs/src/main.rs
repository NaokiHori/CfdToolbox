use komainu::{Document, Expr};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let template = std::fs::read_to_string("template/README.komainu.md")?;
    let mut doc = Document::new(template);
    doc.function("heaviside", "H", 1);
    doc.operator("vector", 2, |args| Expr::subscript(&args[0], &args[1]));
    doc.operator("step_at", 2, |args| Expr::superscript(&args[0], &args[1]));
    doc.operator("partial", 2, |args| {
        Expr::partial_derivative(&args[0], &args[1])
    });
    doc.operator("fraction", 2, |args| Expr::fraction(&args[0], &args[1]));
    doc.operator("particle_volume_integral", 3, |args| {
        Expr::integral_over(&args[0], &args[1], &args[2])
    });
    doc.symbol("pi", r"\pi");
    doc.symbol("i", r"i");
    doc.symbol("j", r"j");
    doc.symbol("k", r"k");
    doc.symbol("ij", r"ij");
    doc.symbol("levi_civita", r"\epsilon_{ijk}");
    doc.symbol("position", r"x");
    doc.symbol("fluid_density", r"\rho^f");
    doc.symbol("external_acceleration", r"g");
    doc.symbol("cauchy_stress_tensor", r"\sigma");
    doc.symbol("current", r"n");
    doc.symbol("predicted", r"*");
    doc.symbol("second_predicted", r"**");
    doc.symbol("next", r"n + 1");
    doc.symbol("fluid_velocity", r"u");
    doc.symbol("particle_center_position", r"\Xi");
    doc.symbol("particle_center_velocity", r"\Upsilon");
    doc.symbol("particle_center_angular_velocity", r"\Omega");
    doc.symbol("particle_surface_position", r"X");
    doc.symbol("particle_surface_velocity", r"U");
    doc.symbol("time_step_size", r"\Delta t");
    doc.symbol("immersed_object_acceleration", "a");
    doc.symbol("particle_density", r"\rho^p");
    doc.symbol("particle_radius", r"R");
    doc.symbol("particle_mass", r"m^p");
    doc.symbol("particle_moment_of_inertia", r"I^p");
    doc.symbol("particle_volume", r"V^p");
    doc.symbol("particle_volume_fraction", r"\phi");
    let markdown = doc.render()?;
    std::fs::write("README.md", markdown)?;
    Ok(())
}
