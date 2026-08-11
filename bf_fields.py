"""Data-driven registry for the Brain-Folding Simulator parameter form.

Every parameter the old ``BFSimulator.py`` exposed is described here once, as a
``Field``: its label, section, the exact ``Parameters.prm`` key it writes, its
validation rule, the info image/text/link shown when it gets focus, and its
2D/3D preset value. ``BFSimulator.py`` builds the whole UI by iterating this
list, which replaces the ~800 lines of copy-pasted widget / ``check_*`` /
``*_info`` code in the original.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional


# --- sections (render order) ------------------------------------------------
GEOMETRY = "Geometry Parameters"
DIFFUSION = "Advection diffusion Parameters"
STIFFNESS = "Mechanical properties Parameters"
MESH = "Discretization Parameters"
SOLVER = "Numerical solver Parameters"
GROWTH = "Growth Parameters"

SECTION_ORDER = [GEOMETRY, DIFFUSION, STIFFNESS, MESH, SOLVER, GROWTH]

# OSVZ distribution combobox options.
OSVZ_OPTIONS = ["Constant", "Linear-gradient", "Quadratic-gradient", "Random1", "Random2"]


@dataclass
class Info:
    """What the info panel shows when a field (or a specific option) is focused."""
    text: str = ""
    img: Optional[tuple[str, str]] = None      # (light, dark) filenames in Images/
    curve: Optional[tuple[str, str]] = None    # secondary image (OSVZ curves)
    url: Optional[str] = None
    url_text: Optional[str] = None


@dataclass
class Field:
    key: str
    label: str
    section: str
    kind: str = "entry"                         # entry | combo | radio
    prm: Optional[str] = None                   # substring key in Parameters.prm (None = not written)
    options: Optional[list[tuple[str, str]]] = None   # (display, value) for radio/combo
    info: object = None                         # Info, or dict[value -> Info] for per-option widgets
    default_2d: str = ""
    default_3d: str = ""
    required: bool = True                       # part of the "empty field" guard on Run
    validate: Optional[Callable[[str, dict], bool]] = None  # returns True when value is OK
    # presentation metadata for the Parameter Guide / unit column
    unit: str = "—"
    symbol: str = "—"
    hint: str = "—"                             # recommended range, human-readable
    error_msg: str = ""                         # short inline validation message


# --- validation helpers -----------------------------------------------------
def _f(s: str):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _v_vz(v, ctx):
    x = _f(v)
    return x is None or (0.2 <= x <= 0.4)


def _v_svz(v, ctx):
    x = _f(v)
    if x is None:
        return True
    lo = _f(ctx.get("vz_raduis", ""))
    lo = lo if lo is not None else 0.0
    return lo <= x <= 0.5


def _v_cr(v, ctx):
    x = _f(v)
    return x is None or (0.01 <= x <= 0.35)


def _v_mst(v, ctx):
    x = _f(v)
    return x is None or x <= 0.1


def _v_poisson(v, ctx):
    x = _f(v)
    return x is None or (0.0 <= x <= 0.5)


def _v_cmax(v, ctx):
    x = _f(v)
    return x is None or x > 200


def _v_stability(v, ctx):
    x = _f(v)
    return x is None or x <= 0.1


def _v_ck(v, ctx):
    x = _f(v)
    return x is None or x <= 1


def _v_newton(v, ctx):
    x = _f(v)
    return x is None or x >= 3


def _v_tol_u(v, ctx):
    x = _f(v)
    return x is None or x <= 1.0e-4


def _v_tol_c(v, ctx):
    x = _f(v)
    return x is None or x <= 1.0e-4


def _v_tol_update(v, ctx):
    x = _f(v)
    return x is None or x <= 1.0e-3


def _v_delt_t(v, ctx):
    x = _f(v)
    return x is None or (0.0 < x < 1.0)


def _v_pos_int(v, ctx):
    if v == "":
        return True
    try:
        return int(v) != 0
    except ValueError:
        return False


def _v_int(v, ctx):
    if v == "":
        return True
    try:
        int(v)
        return True
    except ValueError:
        return False


def _v_k_growth(v, ctx):
    x = _f(v)
    if x is None:
        return True
    case = ctx.get("case", "")
    if case == "2":
        return x <= 1.0e-3
    if case == "3":
        return x <= 1.0e-4
    return True


# --- info-panel content reused across fields --------------------------------
_MORE = "For more details click here"

FIELDS: list[Field] = [
    # ---------------- Geometry ----------------
    Field("vz_raduis", "Ventricular zone raduis:", GEOMETRY, prm="set Ventricular zone raduis",
          default_2d="0.25", default_3d="0.25", validate=_v_vz,
          info=Info("Ventricular zone radius as a ratio to the initial radius.\nShould take a value between 0.2 and 0.4.",
                    img=("gemotry_vz_light.png", "gemotry_vz_dark.png"))),
    Field("svz_raduis", "Subventricular zone raduis:", GEOMETRY, prm="set Subventricular zone raduis",
          default_2d="0.4", default_3d="0.4", validate=_v_svz,
          info=Info("Inner subventricular zone radius as a ratio to the initial radius.\nShould take a value between the ventricular zone radius and 0.5.",
                    img=("gemotry_isvz_light.png", "gemotry_isvz_dark.png"))),
    Field("cr_thickness", "Cortex thickness:", GEOMETRY, prm="set Cortex thickness",
          default_2d="0.1", default_3d="0.1", validate=_v_cr,
          info=Info("Initial cortex thickness as a ratio to the initial radius.\nShould take a value between 0.01 and 0.35.",
                    img=("gemotry_tc_light.png", "gemotry_tc_dark.png"))),
    Field("intial_raduis", "Initial brain radius:", GEOMETRY, prm="set Initial radius",
          default_2d="2", default_3d="2",
          info=Info("Initial fetal brain radius at gestational week 24 [mm].",
                    img=("gemotry_R_light.png", "gemotry_R_dark.png"))),
    Field("MST_factor", "Mitotic translocation factor:", GEOMETRY, prm="set Mitotic somal translocation factor",
          default_2d="0.02", default_3d="0.02", validate=_v_mst,
          info=Info("Mitotic somal translocation factor of ORGCs.\nShould take a value smaller than 0.1.",
                    img=("gemotry_mst_light.png", "gemotry_mst_dark.png"))),

    # ---------------- Advection-diffusion ----------------
    Field("ridial_rate", "Cell division rate of RGCs:", DIFFUSION, prm="set Cell dvision rate of RGCs",
          default_2d="60", default_3d="600",
          info=Info("Division rate in the ventricular zone [1/(mm²·wk)].\nThis factor mimics RGC proliferation.",
                    img=("adv_dif_eq_G_vz_light.png", "adv_dif_eq_G_vz_dark.png"),
                    url="https://en.wikipedia.org/wiki/Radial_glial_cell", url_text=_MORE)),
    Field("outer_ridial_rate", "Cell division rate of Outer RGCs:", DIFFUSION, prm="set Cell dvision rate of Outer RGCs",
          default_2d="10", default_3d="100",
          info=Info("Division rate in the outer subventricular zone [1/(mm²·wk)].\nThis factor mimics ORGC proliferation.",
                    img=("adv_dif_eq_G_osvz_light.png", "adv_dif_eq_G_osvz_dark.png"))),
    Field("ORG_variation_case", "Distribution of OSVZ proliferation:", DIFFUSION,
          kind="combo", prm="set The OSVZ regional variation",
          options=[(o, o) for o in OSVZ_OPTIONS],
          default_2d="Constant", default_3d="Constant",
          info={
              "Constant": Info("This variable controls the regional variation of ORGCs.\nConstant means no regional variation.",
                               img=("OSVZ_Constant.png", "OSVZ_Constant.png")),
              "Linear-gradient": Info("This variable controls the regional variation of ORGCs.\nThe ORGC division rate increases linearly between angles 0° and 90°.",
                                      img=("OSVZ_Linear_gradient.png", "OSVZ_Linear_gradient.png"),
                                      curve=("Linear_gradient_curve.png", "Linear_gradient_curve_dark.png")),
              "Quadratic-gradient": Info("This variable controls the regional variation of ORGCs.\nThe ORGC division rate increases quadratically between angles 0° and 90°.",
                                         img=("OSVZ_Quadratic_gradient.png", "OSVZ_Quadratic_gradient.png"),
                                         curve=("Quadratic_gradient_curve.png", "Quadratic_gradient_curve_dark.png")),
              "Random1": Info("This variable controls the regional variation of ORGCs.\nThe ORGC division rate varies randomly, as shown in the figure.",
                              img=("OSVZ_Random1.png", "OSVZ_Random1.png"),
                              curve=("OSVZ_Random1_curve.png", "OSVZ_Random1_curve_dark.png")),
              "Random2": Info("This variable controls the regional variation of ORGCs.\nThe ORGC division rate varies randomly, as shown in the figure.",
                              img=("OSVZ_Random2.png", "OSVZ_Random2.png"),
                              curve=("OSVZ_Random2_curve.png", "OSVZ_Random2_curve_dark.png")),
          }),
    Field("intial_dvision", "Initial cell density value:", DIFFUSION, prm="set Cell dvision intial value",
          default_2d="0", default_3d="0",
          info=Info("Initial cell density value in the ventricular zone [1/mm²].")),
    Field("migration_speed", "Cell migration speed:", DIFFUSION, prm="set Cell migration speed",
          default_2d="5", default_3d="5",
          info=Info("Cell migration speed [mm/wk].\nThe cells migrate along RGC fibers, i.e. in the radial direction.",
                    img=("adv_dif_eq_v_light.png", "adv_dif_eq_v_dark.png"))),
    Field("diffusivity", "Diffusivity:", DIFFUSION, prm="set Diffusivity",
          default_2d="0.11", default_3d="0.11",
          info=Info("Diffusivity in the cortex [mm²/wk].\nThis model assumes isotropic diffusion, meaning the\ndiffusion is equal in all directions.",
                    img=("adv_dif_eq_d_light.png", "adv_dif_eq_d_dark.png"))),
    Field("migration_threshold", "Cell migration threshold:", DIFFUSION, prm="set Cell migration threshold",
          default_2d="500", default_3d="500",
          info=Info("Cell migration threshold [1/mm³].",
                    img=("adv_dif_eq_c0_light.png", "adv_dif_eq_c0_dark.png"))),
    Field("HV_exp", "Heaviside function exponent:", DIFFUSION, prm="set Heaviside function exponent",
          default_2d="0.008", default_3d="0.008",
          info=Info("Heaviside function exponent.\nFor a smooth solution, it should take a value smaller than 0.01.",
                    img=("adv_dif_eq_gamma_light.png", "adv_dif_eq_gamma_dark.png"))),

    # ---------------- Mechanical properties ----------------
    Field("stiffness_case", "Cortical stiffness case:", STIFFNESS,
          kind="radio", prm="set The state of the stiffness",
          options=[("Varying", "Varying"), ("Constant", "Constant")],
          default_2d="Varying", default_3d="Varying",
          info={
              "Varying": Info("The state of cortical stiffness. 'Varying' means it has a positive relationship with the cell density.",
                              img=("varying.png", "varying_dark.png")),
              "Constant": Info("The cortical stiffness is held constant."),
          }),
    Field("shear_modulus", "Cortical shear modulus:", STIFFNESS, prm="set The shear modulus of conrtex",
          default_2d="2.07", default_3d="2.07",
          info=Info("The shear modulus of the cortical layer [kPa]. The recommended value from the literature is 2.07 kPa.",
                    img=("strain_Energy_mu_light.png", "strain_Energy_mu_dark.png"))),
    Field("stiffness_ratio", "Stiffness ratio:", STIFFNESS, prm="set The ratio of stiffness",
          default_2d="3", default_3d="3",
          info=Info("The ratio of stiffness between cortex and subcortex. The recommended values are 3 and 5.",
                    img=("strain_Energy_ratio_light.png", "strain_Energy_ratio_dark.png"))),
    Field("poisson_ratio", "Poisson's ratio:", STIFFNESS, prm="set Poisson's ratio",
          default_2d="0.38", default_3d="0.38", validate=_v_poisson,
          info=Info("The Poisson's ratio; it must take a value between 0.0 and 0.5.",
                    img=("strain_Energy_nu_light.png", "strain_Energy_nu_dark.png"),
                    url="https://en.wikipedia.org/wiki/Poisson%27s_ratio", url_text=_MORE)),
    Field("max_density", "Maximum cell density:", STIFFNESS, prm="set The max cell density",
          default_2d="700", default_3d="700", validate=_v_cmax,
          info=Info("The maximum cell density. Above this value the stiffness becomes constant and equal to the\nvalue set in the cortical shear modulus. Here c_min is set to 200, so c_max must be greater than 200.",
                    img=("varying_cmax_light.png", "varying_cmax_dark.png"))),

    # ---------------- Discretization ----------------
    Field("case", "Geometry:", MESH, kind="radio",
          options=[("2D", "2"), ("3D", "3")],
          default_2d="2", default_3d="3", required=False,
          info={
              "2": Info("", img=("2D_light.png", "2D_dark.png")),
              "3": Info("", img=("3D_light.png", "3D_dark.png")),
          }),
    Field("refinement", "Number global refinements:", MESH, prm="set Number global refinements",
          default_2d="3", default_3d="2", validate=_v_int,
          info=Info("The number of global mesh refinements. Increasing this value makes the mesh\nfiner and the solution smoother, but the solving time longer.\nThe recommended value is 3 for the 2D case and 2 for the 3D case.",
                    img=("ref_2d_light.png", "ref_2d_dark.png"))),
    Field("degree", "Polynomial degree:", MESH, prm="set Poly degree",
          default_2d="2", default_3d="2", validate=_v_pos_int,
          info=Info("The polynomial degree of the finite-element shape functions, which approximate\nthe solution. Increasing the degree makes the shape functions smoother and the\nsolution more accurate, but increases the solving time.",
                    url="https://en.wikipedia.org/wiki/Hp-FEM", url_text=_MORE)),
    Field("total_time", "Total simulation time:", MESH, prm="set Total time",
          default_2d="1000", default_3d="1000", validate=_v_pos_int,
          info=Info("Total run time. If you do not know the exact time, enter 1000; the solver\nwill then stop automatically when it reaches the mechanical instability point.")),
    Field("delt_t", "Time step size:", MESH, prm="set Time step size",
          default_2d="0.1", default_3d="0.1", validate=_v_delt_t,
          info=Info("Time-step size; it should take a value smaller than 1.0.")),
    Field("stability_con", "Stabilization constant:", MESH, prm="set Stabilization constant",
          default_2d="0.0335", default_3d="0.0335", validate=_v_stability,
          info=Info("Stabilization constant β of the advection–diffusion equation.\nThis model applies a numerical stabilization method.\nThis value should be smaller than 0.1. The recommended value is 0.03334.",
                    url="https://www.dealii.org/current/doxygen/deal.II/step_31.html", url_text=_MORE)),
    Field("c_k", "c_k factor:", MESH, prm="set c_k factor",
          default_2d="0.33334", default_3d="0.33334", validate=_v_ck,
          info=Info("c_k factor that satisfies the CFL condition. This value matters only when the Newton–Raphson\nmethod fails to converge. In that case, the solver re-solves the unconverged time step\nusing a smaller time-step size determined by c_k.",
                    url="https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition", url_text=_MORE)),

    # ---------------- Numerical solver ----------------
    Field("nonlinear_it", "Maximum Newton iterations:", SOLVER, prm="set Max number newton iterations",
          default_2d="8", default_3d="8", required=False, validate=_v_newton,
          info=Info("The maximum number of nonlinear iterations allowed.\nThe Newton–Raphson method is used to solve the nonlinear problem.",
                    url="https://en.wikipedia.org/wiki/Newton%27s_method", url_text=_MORE)),
    Field("tol_u", "Tolerance residual deformation:", SOLVER, prm="set Tolerance residual deformation",
          default_2d="1.0e-8", default_3d="1.0e-8", required=False, validate=_v_tol_u,
          info=Info("Force-residual error tolerance. The recommended value is 1.0e-8.\nThe largest allowed value is 1.0e-4.")),
    Field("tol_c", "Tolerance residual cell density:", SOLVER, prm="set Tolerance residual diffusion",
          default_2d="1.0e-8", default_3d="1.0e-8", required=False, validate=_v_tol_c,
          info=Info("Advection–diffusion residual error tolerance. The recommended value is 1.0e-8.\nThe largest allowed value is 1.0e-4.")),
    Field("update_u", "Tolerance update:", SOLVER, prm="set Tolerance update",
          default_2d="1.0e-4", default_3d="1.0e-4", required=False, validate=_v_tol_update,
          info=Info("Displacement and cell-density update error tolerance. The recommended value is 1.0e-4.\nThe largest allowed value is 1.0e-3.")),
    Field("solver_type", "Linear solver type:", SOLVER,
          kind="radio", prm="set Linear solver type",
          options=[("Direct", "Direct"), ("CG", "CG")],
          default_2d="Direct", default_3d="Direct", required=False,
          info={
              "Direct": Info("The type of solver used for the linear system.\nIf you choose the CG solver, you must enter the number of linear-solver iterations.",
                             url="https://en.wikipedia.org/wiki/Conjugate_gradient_method", url_text=_MORE),
              "CG": Info("The type of solver used for the linear system.\nIf you choose the CG solver, you must enter the number of linear-solver iterations.",
                         url="https://en.wikipedia.org/wiki/Conjugate_gradient_method", url_text=_MORE),
          }),
    Field("linear_it", "Iterations linear solver:", SOLVER, prm="set Multiplier max iterations linear solver",
          default_2d="", default_3d="", required=False,
          info=Info("The maximum number of iterations for the CG linear solver.")),

    # ---------------- Growth ----------------
    Field("k_growth", "Growth rate:", GROWTH, prm="set Growth rate",
          default_2d="4.7e-4", default_3d="4.7e-5", validate=_v_k_growth,
          info=Info("Growth-rate factor. This constant coefficient controls the amount of isotropic growth in\nthe subcortical layer. For numerical stability, it must be smaller than 1.0e-3 for the\n2D case and 1.0e-4 for the 3D case.",
                    img=("growth_eqautions_ks_light.png", "growth_eqautions_ks_dark.png"))),
    Field("growth_ratio", "Growth ratio:", GROWTH, prm="set  Growth ratio",
          default_2d="1.5", default_3d="2",
          info=Info("Growth ratio. This ratio controls the amount of tangential versus radial growth in the\ncortical layer. Increasing it increases the variation between tangential and radial growth.\nThe recommended values are 1.5 and 3.",
                    img=("growth_eqautions_ratio_light.png", "growth_eqautions_ratio_dark.png"))),
    Field("growth_exp", "Growth exponent:", GROWTH, prm="set Growth exponent",
          default_2d="1.65", default_3d="1.65",
          info=Info("Growth exponent.",
                    img=("growth_eqautions_exp_light.png", "growth_eqautions_exp_dark.png"))),
]

FIELDS_BY_KEY = {f.key: f for f in FIELDS}

# --- presentation metadata: unit, symbol, recommended range, inline error -----
# (unit, symbol, hint, error_msg) keyed by field. Symbols use plain text.
META: dict[str, tuple[str, str, str, str]] = {
    "vz_raduis":          ("ratio", r"r_{VZ}", "0.2 – 0.4", "Must be between 0.2 and 0.4."),
    "svz_raduis":         ("ratio", r"r_{ISVZ}", "≥ r_vz, ≤ 0.5", "Must be ≥ ventricular radius and ≤ 0.5."),
    "cr_thickness":       ("ratio", r"t_c", "0.01 – 0.35", "Must be between 0.01 and 0.35."),
    "intial_raduis":      ("mm", "R", "> 0", ""),
    "MST_factor":         ("ratio", r"m_{mst}", "≤ 0.1", "Should be smaller than 0.1."),
    "ridial_rate":        ("1/(mm²·wk)", r"G_{VZ}", "—", ""),
    "outer_ridial_rate":  ("1/(mm²·wk)", r"G_{OSVZ}", "—", ""),
    "ORG_variation_case": ("—", "—", "5 modes", "Choose a distribution."),
    "intial_dvision":     ("1/mm²", "—", "—", ""),
    "migration_speed":    ("mm/wk", "v", "—", ""),
    "diffusivity":        ("mm²/wk", r"d^{cc}", "—", ""),
    "migration_threshold": ("1/mm³", r"c_0", "—", ""),
    "HV_exp":             ("—", r"\gamma", "< 0.01", "Should be smaller than 0.01."),
    "stiffness_case":     ("—", "—", "Varying / Constant", ""),
    "shear_modulus":      ("kPa", r"\mu", "≈ 2.07", ""),
    "stiffness_ratio":    ("—", r"\beta_{\mu}", "3 – 5", ""),
    "poisson_ratio":      ("—", r"\nu", "0.0 – 0.5", "Must be between 0.0 and 0.5."),
    "max_density":        ("1/mm³", r"c_{max}", "> 200", "Must be greater than 200."),
    "case":               ("—", "—", "2D / 3D", ""),
    "refinement":         ("—", "—", "2D:3 · 3D:2", "Must be an integer."),
    "degree":             ("—", "—", "≥ 1", "Must be a non-zero integer."),
    "total_time":         ("s", "—", "> 0", "Must be a non-zero integer."),
    "delt_t":             ("s", r"\Delta t", "< 1.0", "Must be between 0 and 1.0."),
    "stability_con":      ("—", r"\beta", "< 0.1", "Should be smaller than 0.1."),
    "c_k":                ("—", r"c_k", "≤ 1", "Must be ≤ 1."),
    "nonlinear_it":       ("—", "—", "≥ 3", "Must be at least 3."),
    "tol_u":              ("—", "—", "≤ 1e-4", "Largest allowed value is 1e-4."),
    "tol_c":              ("—", "—", "≤ 1e-4", "Largest allowed value is 1e-4."),
    "update_u":           ("—", "—", "≤ 1e-3", "Largest allowed value is 1e-3."),
    "solver_type":        ("—", "—", "Direct / CG", ""),
    "linear_it":          ("—", "—", "CG only", ""),
    "k_growth":           ("—", r"k_s", "2D ≤1e-3 · 3D ≤1e-4", "Exceeds the stability limit for this case."),
    "growth_ratio":       ("—", r"\beta_k", "1.5 – 3", ""),
    "growth_exp":         ("—", r"\alpha", "—", ""),
}
for _fld in FIELDS:
    _u, _s, _h, _e = META.get(_fld.key, ("—", "—", "—", "Invalid value."))
    _fld.unit, _fld.symbol, _fld.hint, _fld.error_msg = _u, _s, _h, _e

# --- per-category identity: accent key + line-icon SVG path ------------------
# accent keys resolve to hex in bf_style.ACCENTS. Icons are 24x24 stroke paths.
SECTION_ACCENT: dict[str, str] = {
    GEOMETRY: "geo", DIFFUSION: "adv", STIFFNESS: "mech",
    MESH: "disc", SOLVER: "solv", GROWTH: "grow",
}
SECTION_ICON: dict[str, str] = {
    # Cube
    GEOMETRY:  '<path d="M12 2l9 5v10l-9 5-9-5V7z"/><path d="M12 12l9-5M12 12v10M12 12L3 7"/>',
    # Flowing lines
    DIFFUSION: '<path d="M3 8c3 0 3 3 6 3s3-3 6-3 3 3 6 3M3 15c3 0 3 3 6 3s3-3 6-3 3 3 6 3"/>',
    # Spring / coil
    STIFFNESS: '<path d="M2 12h2M20 12h2"/><path d="M4 12c0-4 2.7-4 2.7 0s2.7 4 2.7 0 2.7-4 2.7 0 2.7 4 2.7 0"/>',
    # Mesh / grid
    MESH:      '<rect x="4" y="4" width="16" height="16" rx="1"/><path d="M4 10h16M4 15h16M10 4v16M15 4v16"/>',
    # Sigma (Σ)
    SOLVER:    '<path d="M17 5H7l6 7-6 7h10"/>',
    # Rising arrow / chart
    GROWTH:    '<path d="M3 17l6-6 4 4 8-8"/><path d="M21 11V7h-4"/>',
}

# Open-book icon for the Parameter Guide panel header.
GUIDE_ICON = ('<path d="M12 6C10 4.5 7 4.5 4 6v12c3-1.5 6-1.5 8 0 '
              '2-1.5 5-1.5 8 0V6c-3-1.5-6-1.5-8 0z"/><path d="M12 6v12"/>')
SECTION_SHORT: dict[str, str] = {
    GEOMETRY: "Geometry", DIFFUSION: "Advection–Diffusion", STIFFNESS: "Mechanical Properties",
    MESH: "Discretization", SOLVER: "Numerical Solver", GROWTH: "Growth",
}


# --- static info panels reachable from the buttons --------------------------
ABOUT_AUTHOR = Info(
    "Mohammad Saeed Zarzor\n\n"
    "- PhD candidate in the field of Biomechanics.\n"
    "- Scientific employee in Institute of Applied Mechanics,\n"
    "  Friedrich-Alexander-University Erlangen-Nürnberg.\n"
    "- Master of Science in Computational Engineering from FAU University.\n"
    "- Bachelor of Mechanical Engineering from Damascus University.",
    img=("my_photo.png", "my_photo.png"),
    url="https://www.ltm.tf.fau.eu/person/zarzor-mohammad-saeed-m-sc/", url_text="contact details")

ABOUT_PROGRAM = Info(
    "This work is part of the BRAINIACS Project at the Institute of Applied\n"
    "Mechanics, Friedrich-Alexander-University Erlangen-Nürnberg, under the\n"
    "supervision of Dr. Silvia Budday and in cooperation with Prof. Dr. med\n"
    "Ingmar Blümcke from Neuropathological Institute, University Hospitals\n"
    "Erlangen. We gratefully acknowledge the funding by the Deutsche\n"
    "Forschungsgemeinschaft (DFG, German Research Foundation).\n"
    "This work is based on the following paper:",
    img=("Logo_BRAINIACS.png", "Logo_BRAINIACS.png"),
    url="https://www.biorxiv.org/content/10.1101/2022.09.25.509401v1.abstract",
    url_text="Exploring the role of the outer subventricular zone during cortical folding through a physics-based model")

COPYRIGHT = Info(
    "The copyright holder for this preprint is the author/funder, who has granted bioRxiv a license\n"
    "to display the preprint in perpetuity. It is made available under a Copyright:",
    url="https://creativecommons.org/licenses/by/4.0/",
    url_text="CC-BY 4.0 International license.")

BACKGROUND = ("Untitled-3.png", "Untitled-3-dark.png")
