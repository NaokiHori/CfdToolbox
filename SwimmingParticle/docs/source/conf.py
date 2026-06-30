project = "Swimming Particle"
author = "Naoki Hori"
copyright = f"2025, {author}"

html_theme = "alabaster"
html_static_path = ["_static"]
html_theme_options = {
    "fixed_sidebar": "false",
    "navigation_with_keys": "true",
    "nosidebar": "false",
}

macros = dict()
macros["vr"] = ["r", 0]
macros["vt"] = ["\\theta", 0]
macros["ur"] = ["u_{\\vr}", 0]
macros["ut"] = ["u_{\\vt}", 0]
macros["nr"] = ["N_{\\vr}", 0]
macros["nt"] = ["N_{\\vt}", 0]
macros["dim"] = ["\\tilde{#1}", 1]
macros["vat"] = ["\\left. #1 \\right|_{#2}", 2]
macros["der"] = ["\\frac{d^{#1} {#2}}{d {#3}^{#1}}", 3]
macros["pder"] = ["\\frac{\\partial^{#1} {#2}}{\\partial {#3}^{#1}}", 3]
macros["dder"] = ["\\frac{\\delta^{#1} {#2}}{\\delta {#3}^{#1}}", 3]
macros["expp"] = ["\\exp \\left( {#1} \\right)", 1]

# discrete advective terms in scalar transport
macros["ave"] = ["\\overline{#1}^{#2}", 2]
macros["dif"] = ["{\\delta_{#2} {#1}}", 2]
macros["gcs"] = ["{\\xi^{#1}}", 1]
macros["sfact"] = ["h_{\\gcs{#1}}", 1]
macros["vel"] = ["{u_{#1}}", 1]
macros["scalar"] = ["c", 0]
macros["dscalaradv"] = [
        "\\frac{1}{J}"
        "\\ave{"
        "  \\frac{J}{\\sfact{#1}}"
        "  \\vel{#1}"
        "  \\dif{#2}{\\gcs{#1}}"
        "}{\\gcs{#1}}"
        , 2
]
# discrete diffusive terms in scalar transport
macros["dscalardif"] = [
        "\\frac{1}{J}"
        "\\dif{"
        "}{\\gcs{#1}}"
        "\\left("
        "  \\frac{J}{\\sfact{#1}}"
        "  \\frac{1}{\\sfact{#1}}"
        "  \\dif{#2}{\\gcs{#1}}"
        "\\right)"
        , 2
]

mathjax3_config = {"tex": {"macros": macros}}
