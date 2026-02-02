"""
bench.py — Bench chemistry calculators and reference databases.

Pure-computation tools for lab work: molarity, dilution, yield, and mass calculations,
plus comprehensive reference data for named reactions, protecting groups, workup
procedures, solvents, cooling baths, TLC stains, and column chromatography.

No API keys required. Works out of the box.
"""

from __future__ import annotations

import math
import re
from typing import Optional

# =============================================================================
# Constants — Molecular weight helpers
# =============================================================================

ATOMIC_WEIGHTS = {
    "H": 1.008, "He": 4.003, "Li": 6.941, "Be": 9.012, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.086, "P": 30.974,
    "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.64, "As": 74.922, "Se": 78.96, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.96, "Ru": 101.07, "Rh": 102.906, "Pd": 106.42,
    "Ag": 107.868, "Cd": 112.411, "In": 114.818, "Sn": 118.710, "Sb": 121.760,
    "Te": 127.60, "I": 126.904, "Xe": 131.293, "Cs": 132.905, "Ba": 137.327,
    "La": 138.905, "Ce": 140.116, "Pr": 140.908, "Nd": 144.242, "Sm": 150.36,
    "Eu": 151.964, "Gd": 157.25, "Tb": 158.925, "Dy": 162.500, "Ho": 164.930,
    "Er": 167.259, "Tm": 168.934, "Yb": 173.054, "Lu": 174.967, "Hf": 178.49,
    "Ta": 180.948, "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217,
    "Pt": 195.084, "Au": 196.967, "Hg": 200.59, "Tl": 204.383, "Pb": 207.2,
    "Bi": 208.980, "Th": 232.038, "U": 238.029,
}


def parse_formula(formula: str) -> float:
    """Parse a molecular formula string and return molecular weight.
    Handles nested parentheses, e.g. Ca(OH)2, Mg3(PO4)2."""
    formula = formula.strip()
    # Recursive approach: handle parentheses first
    while "(" in formula:
        # Find innermost parentheses
        m = re.search(r"\(([^()]+)\)(\d*)", formula)
        if not m:
            break
        inner = m.group(1)
        mult = int(m.group(2)) if m.group(2) else 1
        # Expand: multiply each element count
        inner_mw = _sum_elements(inner)
        # Replace with placeholder
        formula = formula[:m.start()] + f"[{inner_mw * mult}]" + formula[m.end():]

    # Now parse remaining, including placeholders
    total = 0.0
    # Handle bracketed MW values from parenthesis expansion
    for bm in re.finditer(r"\[([\d.]+)\]", formula):
        total += float(bm.group(1))
    formula_clean = re.sub(r"\[[\d.]+\]", "", formula)
    total += _sum_elements(formula_clean)
    return total


def _sum_elements(s: str) -> float:
    """Sum molecular weights of elements in a simple formula (no parentheses)."""
    total = 0.0
    for m in re.finditer(r"([A-Z][a-z]?)(\d*)", s):
        elem = m.group(1)
        count = int(m.group(2)) if m.group(2) else 1
        if elem in ATOMIC_WEIGHTS:
            total += ATOMIC_WEIGHTS[elem] * count
    return total


# =============================================================================
# Calculator: Molarity & Solutions
# =============================================================================

def calc_molarity(
    *,
    concentration: float | None = None,
    mass_grams: float | None = None,
    moles: float | None = None,
    volume_ml: float | None = None,
    volume_l: float | None = None,
    mw: float | None = None,
    formula: str | None = None,
) -> dict:
    """Multi-purpose molarity/solution calculator.

    Provide 2-3 of: concentration (M), mass (g), moles, volume (mL or L), MW.
    Solves for unknowns using M = n/V, n = mass/MW.
    """
    # Resolve volume to liters
    vol_l = None
    if volume_l is not None:
        vol_l = volume_l
    elif volume_ml is not None:
        vol_l = volume_ml / 1000.0

    # Resolve MW from formula if needed
    if mw is None and formula:
        mw = parse_formula(formula)

    # Resolve moles from mass and MW
    if moles is None and mass_grams is not None and mw is not None and mw > 0:
        moles = mass_grams / mw

    # Resolve mass from moles and MW
    if mass_grams is None and moles is not None and mw is not None:
        mass_grams = moles * mw

    # M = n / V
    if concentration is None and moles is not None and vol_l is not None and vol_l > 0:
        concentration = moles / vol_l
    if moles is None and concentration is not None and vol_l is not None:
        moles = concentration * vol_l
    if vol_l is None and concentration is not None and moles is not None and concentration > 0:
        vol_l = moles / concentration

    # Re-derive mass if we now have moles
    if mass_grams is None and moles is not None and mw is not None:
        mass_grams = moles * mw
    if moles is None and mass_grams is not None and mw is not None and mw > 0:
        moles = mass_grams / mw

    return {
        "molarity_M": _r(concentration),
        "moles_mol": _r(moles),
        "mass_g": _r(mass_grams),
        "volume_mL": _r(vol_l * 1000) if vol_l is not None else None,
        "volume_L": _r(vol_l),
        "mw_g_per_mol": _r(mw),
    }


def _r(v, dp=6):
    if v is None:
        return None
    return round(v, dp)


# =============================================================================
# Calculator: Dilution (C₁V₁ = C₂V₂)
# =============================================================================

def calc_dilution(
    *,
    c1: float | None = None,
    v1: float | None = None,
    c2: float | None = None,
    v2: float | None = None,
    c1_unit: str = "M",
    v1_unit: str = "mL",
    c2_unit: str = "M",
    v2_unit: str = "mL",
) -> dict:
    """C₁V₁ = C₂V₂ dilution calculator. Solves for the one missing value."""
    # Normalize concentrations to M
    c1_m = _to_molar(c1, c1_unit) if c1 is not None else None
    c2_m = _to_molar(c2, c2_unit) if c2 is not None else None
    v1_ml = _to_ml(v1, v1_unit) if v1 is not None else None
    v2_ml = _to_ml(v2, v2_unit) if v2 is not None else None

    nones = sum(1 for x in [c1_m, v1_ml, c2_m, v2_ml] if x is None)
    if nones != 1:
        return {"error": f"Provide exactly 3 of 4 values (c1, v1, c2, v2). Got {4 - nones} values."}

    if c1_m is None:
        c1_m = (c2_m * v2_ml) / v1_ml
    elif v1_ml is None:
        v1_ml = (c2_m * v2_ml) / c1_m
    elif c2_m is None:
        c2_m = (c1_m * v1_ml) / v2_ml
    elif v2_ml is None:
        v2_ml = (c1_m * v1_ml) / c2_m

    solvent_to_add = v2_ml - v1_ml if v1_ml is not None and v2_ml is not None else None

    return {
        "c1": _r(c1_m), "c1_unit": "M",
        "v1_mL": _r(v1_ml),
        "c2": _r(c2_m), "c2_unit": "M",
        "v2_mL": _r(v2_ml),
        "solvent_to_add_mL": _r(solvent_to_add) if solvent_to_add and solvent_to_add > 0 else None,
        "dilution_factor": _r(c1_m / c2_m) if c1_m and c2_m else None,
    }


def _to_molar(val: float, unit: str) -> float:
    u = unit.lower().strip()
    if u in ("m", "mol/l"):
        return val
    if u in ("mm", "mmol/l"):
        return val / 1000.0
    if u in ("um", "µm", "μm", "umol/l"):
        return val / 1_000_000.0
    if u in ("nm", "nmol/l"):
        return val / 1_000_000_000.0
    if u in ("mg/ml", "g/l"):
        return val  # user needs MW to convert; treat as direct
    return val  # assume M


def _to_ml(val: float, unit: str) -> float:
    u = unit.lower().strip()
    if u in ("ml",):
        return val
    if u in ("l",):
        return val * 1000.0
    if u in ("ul", "µl", "μl"):
        return val / 1000.0
    return val  # assume mL


# =============================================================================
# Calculator: Reaction Mass / Equivalents
# =============================================================================

def calc_reaction_mass(
    reagents: list[dict],
) -> dict:
    """Calculate masses/volumes needed for a reaction.

    Each reagent dict should have:
      - name: str
      - mw: float (g/mol) OR formula: str
      - equiv: float (equivalents, default 1.0)
      - AND one of: mass_g, moles, volume_ml + density
    Optional per reagent:
      - density: float (g/mL) for liquid reagents
      - molarity: float (M) for solutions
      - volume_ml: float

    The first reagent with mass_g or moles is the reference.
    """
    results = []
    ref_moles = None

    # First pass: resolve MW and find reference moles
    for r in reagents:
        mw = r.get("mw")
        if mw is None and r.get("formula"):
            mw = parse_formula(r["formula"])
        r["_mw"] = mw

        equiv = r.get("equiv", 1.0)
        r["_equiv"] = equiv

        # Find reference compound (first with defined moles)
        if ref_moles is None:
            if r.get("moles"):
                ref_moles = r["moles"] / equiv
            elif r.get("mass_g") and mw and mw > 0:
                ref_moles = (r["mass_g"] / mw) / equiv
            elif r.get("volume_ml") and r.get("density") and mw and mw > 0:
                mass = r["volume_ml"] * r["density"]
                ref_moles = (mass / mw) / equiv
            elif r.get("volume_ml") and r.get("molarity"):
                mol = r["volume_ml"] / 1000.0 * r["molarity"]
                ref_moles = mol / equiv

    if ref_moles is None:
        return {"error": "Could not determine reference moles. Provide mass_g or moles for at least one reagent."}

    # Second pass: compute all quantities
    for r in reagents:
        mw = r["_mw"]
        equiv = r["_equiv"]
        moles_needed = ref_moles * equiv

        entry = {
            "name": r.get("name", "?"),
            "mw": _r(mw),
            "equiv": equiv,
            "moles": _r(moles_needed, 6),
            "mmoles": _r(moles_needed * 1000, 4),
        }

        if mw and mw > 0:
            entry["mass_g"] = _r(moles_needed * mw, 4)
            entry["mass_mg"] = _r(moles_needed * mw * 1000, 2)

        density = r.get("density")
        if density and mw and mw > 0:
            vol = (moles_needed * mw) / density
            entry["volume_mL"] = _r(vol, 4)
            entry["volume_uL"] = _r(vol * 1000, 1)

        molarity = r.get("molarity")
        if molarity and molarity > 0:
            vol_l = moles_needed / molarity
            entry["volume_mL"] = _r(vol_l * 1000, 4)
            entry["volume_uL"] = _r(vol_l * 1_000_000, 1)

        results.append(entry)

    return {
        "ref_moles_1_equiv": _r(ref_moles, 6),
        "reagents": results,
    }


# =============================================================================
# Calculator: Yield
# =============================================================================

def calc_yield(
    *,
    actual_mass_g: float | None = None,
    theoretical_mass_g: float | None = None,
    limiting_moles: float | None = None,
    limiting_mass_g: float | None = None,
    limiting_mw: float | None = None,
    limiting_formula: str | None = None,
    product_mw: float | None = None,
    product_formula: str | None = None,
    stoich_coeff_product: float = 1.0,
    stoich_coeff_limiting: float = 1.0,
    percent_yield: float | None = None,
) -> dict:
    """Calculate theoretical yield and percent yield."""
    # Resolve MWs
    if limiting_mw is None and limiting_formula:
        limiting_mw = parse_formula(limiting_formula)
    if product_mw is None and product_formula:
        product_mw = parse_formula(product_formula)

    # Resolve moles of limiting reagent
    if limiting_moles is None and limiting_mass_g is not None and limiting_mw:
        limiting_moles = limiting_mass_g / limiting_mw

    # Theoretical moles of product
    theo_moles = None
    if limiting_moles is not None:
        theo_moles = limiting_moles * (stoich_coeff_product / stoich_coeff_limiting)

    # Theoretical mass
    if theoretical_mass_g is None and theo_moles is not None and product_mw:
        theoretical_mass_g = theo_moles * product_mw

    # Percent yield
    if percent_yield is None and actual_mass_g is not None and theoretical_mass_g:
        percent_yield = (actual_mass_g / theoretical_mass_g) * 100.0

    # Recover actual mass from percent yield
    if actual_mass_g is None and percent_yield is not None and theoretical_mass_g:
        actual_mass_g = (percent_yield / 100.0) * theoretical_mass_g

    # Atom economy (if we have both MWs)
    atom_economy = None
    if product_mw and limiting_mw:
        atom_economy = (product_mw * stoich_coeff_product) / (limiting_mw * stoich_coeff_limiting) * 100.0
        if atom_economy > 100:
            atom_economy = None  # doesn't make sense w/o all reactants

    return {
        "limiting_moles": _r(limiting_moles),
        "theoretical_moles_product": _r(theo_moles),
        "theoretical_mass_g": _r(theoretical_mass_g, 4),
        "actual_mass_g": _r(actual_mass_g, 4),
        "percent_yield": _r(percent_yield, 2),
        "product_mw": _r(product_mw),
        "limiting_mw": _r(limiting_mw),
    }


# =============================================================================
# Calculator: Concentration Interconversion
# =============================================================================

def calc_concentration(
    *,
    value: float,
    from_unit: str,
    to_unit: str,
    mw: float | None = None,
    formula: str | None = None,
    density_solution: float = 1.0,  # g/mL
) -> dict:
    """Convert between concentration units.

    Supported units: M, mM, uM, nM, mg/mL, g/L, ug/mL, %w/v, %w/w, %v/v, ppm, ppb, mol/L, mmol/L.
    MW is required for conversions between molar and mass-based units.
    """
    if mw is None and formula:
        mw = parse_formula(formula)

    # Normalize everything to mol/L (M) if molar, or mg/mL if mass-based
    # We convert to an intermediate: mg/mL AND mol/L where possible

    u_from = from_unit.lower().strip().replace(" ", "")
    u_to = to_unit.lower().strip().replace(" ", "")

    # Convert from_unit → mg/mL
    mg_per_ml = None
    mol_per_l = None

    if u_from in ("m", "mol/l"):
        mol_per_l = value
        if mw:
            mg_per_ml = value * mw  # M * g/mol = g/L = mg/mL
    elif u_from in ("mm", "mmol/l"):
        mol_per_l = value / 1000.0
        if mw:
            mg_per_ml = mol_per_l * mw
    elif u_from in ("um", "µm", "μm", "umol/l"):
        mol_per_l = value / 1e6
        if mw:
            mg_per_ml = mol_per_l * mw
    elif u_from in ("nm", "nmol/l"):
        mol_per_l = value / 1e9
        if mw:
            mg_per_ml = mol_per_l * mw
    elif u_from in ("mg/ml", "g/l"):
        mg_per_ml = value
        if mw and mw > 0:
            mol_per_l = value / mw
    elif u_from in ("ug/ml", "µg/ml", "mg/l"):
        mg_per_ml = value / 1000.0
        if mw and mw > 0:
            mol_per_l = mg_per_ml / mw
    elif u_from in ("ng/ml", "ug/l", "µg/l"):
        mg_per_ml = value / 1e6
        if mw and mw > 0:
            mol_per_l = mg_per_ml / mw
    elif u_from in ("%w/v", "%(w/v)"):
        mg_per_ml = value * 10.0  # 1% w/v = 1 g/100 mL = 10 mg/mL
        if mw and mw > 0:
            mol_per_l = (mg_per_ml) / mw
    elif u_from in ("%w/w", "%(w/w)"):
        mg_per_ml = value * 10.0 * density_solution
        if mw and mw > 0:
            mol_per_l = mg_per_ml / mw
    elif u_from in ("ppm",):
        mg_per_ml = value / 1000.0  # 1 ppm = 1 mg/L = 0.001 mg/mL
        if mw and mw > 0:
            mol_per_l = mg_per_ml / mw
    elif u_from in ("ppb",):
        mg_per_ml = value / 1e6
        if mw and mw > 0:
            mol_per_l = mg_per_ml / mw
    else:
        return {"error": f"Unknown source unit: {from_unit}"}

    # Convert mg/mL or mol/L → target unit
    result = None

    if u_to in ("m", "mol/l"):
        result = mol_per_l
    elif u_to in ("mm", "mmol/l"):
        result = mol_per_l * 1e3 if mol_per_l is not None else None
    elif u_to in ("um", "µm", "μm", "umol/l"):
        result = mol_per_l * 1e6 if mol_per_l is not None else None
    elif u_to in ("nm", "nmol/l"):
        result = mol_per_l * 1e9 if mol_per_l is not None else None
    elif u_to in ("mg/ml", "g/l"):
        result = mg_per_ml
    elif u_to in ("ug/ml", "µg/ml", "mg/l"):
        result = mg_per_ml * 1e3 if mg_per_ml is not None else None
    elif u_to in ("ng/ml", "ug/l", "µg/l"):
        result = mg_per_ml * 1e6 if mg_per_ml is not None else None
    elif u_to in ("%w/v", "%(w/v)"):
        result = mg_per_ml / 10.0 if mg_per_ml is not None else None
    elif u_to in ("%w/w", "%(w/w)"):
        result = (mg_per_ml / (10.0 * density_solution)) if mg_per_ml is not None else None
    elif u_to in ("ppm",):
        result = mg_per_ml * 1000.0 if mg_per_ml is not None else None
    elif u_to in ("ppb",):
        result = mg_per_ml * 1e6 if mg_per_ml is not None else None
    else:
        return {"error": f"Unknown target unit: {to_unit}"}

    if result is None:
        return {"error": f"Cannot convert from {from_unit} to {to_unit} without molecular weight."}

    return {
        "input_value": value,
        "input_unit": from_unit,
        "result_value": _r(result, 6),
        "result_unit": to_unit,
        "mw_used": _r(mw) if mw else None,
        "density_used": density_solution if u_from in ("%w/w", "%(w/w)") or u_to in ("%w/w", "%(w/w)") else None,
        "intermediate_mg_per_mL": _r(mg_per_ml, 6) if mg_per_ml is not None else None,
        "intermediate_mol_per_L": _r(mol_per_l, 9) if mol_per_l is not None else None,
    }


# =============================================================================
# Named Reactions Database
# =============================================================================

NAMED_REACTIONS: dict[str, dict] = {}

def _nr(name, *, aliases=None, category="", type_="", summary="",
        conditions="", substrate="", limitations="", mechanism="",
        refs="", related="", functional_groups=""):
    """Register a named reaction."""
    key = name.lower()
    NAMED_REACTIONS[key] = {
        "name": name,
        "aliases": aliases or [],
        "category": category,
        "type": type_,
        "summary": summary,
        "conditions": conditions,
        "substrate": substrate,
        "limitations": limitations,
        "mechanism": mechanism,
        "refs": refs,
        "related": related,
        "functional_groups": functional_groups,
    }
    # Also register aliases
    for a in (aliases or []):
        NAMED_REACTIONS[a.lower()] = NAMED_REACTIONS[key]


# --- C-C Bond Formation: Cross-Coupling ---

_nr("Suzuki Coupling",
    aliases=["Suzuki-Miyaura", "Suzuki reaction"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of organoboron compounds (boronic acids/esters) with organic halides or triflates.",
    conditions="Pd(PPh3)4 or Pd(dppf)Cl2 (0.5-5 mol%), base (K2CO3, Cs2CO3, K3PO4), solvent (THF, dioxane, DMF/H2O), 60-100°C. Aqueous conditions possible.",
    substrate="Aryl/vinyl boronic acids + aryl/vinyl halides (I > Br >> Cl) or triflates. sp2-sp2 coupling preferred. sp2-sp3 possible with specialized ligands (SPhos, XPhos).",
    limitations="sp3-sp3 coupling difficult. β-Hydride elimination with alkyl halides. Protodeboronation can compete. Homo-coupling side reaction.",
    mechanism="Oxidative addition → transmetalation (base-assisted) → reductive elimination. Pd(0)/Pd(II) cycle.",
    refs="Miyaura, Suzuki, Chem. Rev. 1995, 95, 2457. Martin, Buchwald, Acc. Chem. Res. 2008, 41, 1461.",
    related="Heck, Negishi, Stille, Kumada, Hiyama",
    functional_groups="aryl halide, boronic acid, boronate ester")

_nr("Heck Reaction",
    aliases=["Mizoroki-Heck", "Heck coupling"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of aryl/vinyl halides with alkenes to form substituted alkenes.",
    conditions="Pd(OAc)2 or PdCl2 (1-5 mol%), ligand (PPh3 or P(o-tol)3), base (Et3N, K2CO3, NaOAc), DMF or DMA, 80-120°C.",
    substrate="Aryl/vinyl halides + alkenes (electron-poor preferred). Regioselectivity: addition at less hindered end of alkene.",
    limitations="Requires β-hydrogens on alkene. Regioselectivity can be problematic. Double bond isomerization possible.",
    mechanism="Oxidative addition → syn-migratory insertion into alkene → β-hydride elimination. Pd(0)/Pd(II) cycle.",
    refs="Heck, Nolley, J. Org. Chem. 1972, 37, 2320. Beletskaya, Cheprakov, Chem. Rev. 2000, 100, 3009.",
    related="Suzuki, Negishi, Stille, Fujiwara-Moritani",
    functional_groups="aryl halide, alkene")

_nr("Negishi Coupling",
    aliases=["Negishi reaction"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd or Ni-catalyzed coupling of organozinc reagents with organic halides. Excellent functional group tolerance.",
    conditions="Pd(PPh3)4 or Pd(dba)2/ligand, THF, 0°C to RT. Organozinc prepared from RLi or RMgX + ZnCl2.",
    substrate="RZnX + R'X (X = I, Br, Cl, OTf). Works with sp, sp2, and sp3 centers. Tolerates esters, nitriles, ketones.",
    limitations="Organozinc reagents are moisture-sensitive (less than Grignard). Requires preparation of organozinc species.",
    mechanism="Oxidative addition → transmetalation → reductive elimination. Pd(0)/Pd(II) cycle.",
    refs="Negishi, King, Okukado, J. Org. Chem. 1977, 42, 1821. Nobel Prize 2010.",
    related="Suzuki, Stille, Kumada",
    functional_groups="organozinc, aryl halide")

_nr("Stille Coupling",
    aliases=["Stille reaction", "Migita-Kosugi-Stille"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of organostannanes (organotin) with organic halides. Very broad scope and excellent FG tolerance.",
    conditions="Pd(PPh3)4 or Pd2(dba)3 + AsPh3 (1-5 mol%), DMF, toluene, or NMP, 50-100°C. CuI or CsF as additives improve rate.",
    substrate="R3SnR' + R''X (X = I, Br, Cl, OTf). Works well with vinyl, aryl, acyl halides. Tolerates almost all functional groups.",
    limitations="Toxicity of organotin compounds. Tin residues difficult to remove. Homo-coupling. Slow transmetalation sometimes.",
    mechanism="Oxidative addition → transmetalation → reductive elimination. Pd(0)/Pd(II) cycle.",
    refs="Stille, Angew. Chem. Int. Ed. 1986, 25, 508. Espinet, Echavarren, Angew. Chem. Int. Ed. 2004, 43, 4704.",
    related="Suzuki, Negishi, Heck",
    functional_groups="organostannane, aryl halide")

_nr("Kumada Coupling",
    aliases=["Kumada-Corriu", "Kumada reaction"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Ni or Pd-catalyzed coupling of Grignard reagents with organic halides. First transition metal-catalyzed cross-coupling.",
    conditions="NiCl2(dppf) or Pd(dppf)Cl2, THF or Et2O, 0°C to RT.",
    substrate="RMgX + R'X. Works best for sp2-sp2 coupling.",
    limitations="Limited FG tolerance (Grignard reacts with ketones, esters, etc.). Cannot use protic solvents.",
    mechanism="Oxidative addition → transmetalation → reductive elimination.",
    refs="Kumada, Tamao, Sumitani, J. Am. Chem. Soc. 1972, 94, 4374. Corriu, Masse, J. Chem. Soc. Chem. Commun. 1972, 144.",
    related="Suzuki, Negishi, Grignard reaction",
    functional_groups="Grignard reagent, aryl halide")

_nr("Sonogashira Coupling",
    aliases=["Sonogashira reaction"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd/Cu-catalyzed coupling of terminal alkynes with aryl/vinyl halides to form aryl/vinyl alkynes.",
    conditions="PdCl2(PPh3)2 (2-5 mol%), CuI (5-10 mol%), amine base (Et3N, iPr2NH, piperidine), THF or DMF, RT-80°C.",
    substrate="Terminal alkynes + aryl/vinyl halides or triflates. I > Br > Cl reactivity.",
    limitations="Glaser homo-coupling of alkyne (in presence of O2). Requires anhydrous/degassed conditions. Cu-free variants available.",
    mechanism="Pd cycle: oxidative addition → transmetalation from Cu-acetylide → reductive elimination. Cu cycle: deprotonation of alkyne → Cu-acetylide formation.",
    refs="Sonogashira, Tohda, Hagihara, Tetrahedron Lett. 1975, 4467. Chinchilla, Nájera, Chem. Rev. 2007, 107, 874.",
    related="Suzuki, Heck, Cadiot-Chodkiewicz, Glaser",
    functional_groups="terminal alkyne, aryl halide")

_nr("Buchwald-Hartwig Amination",
    aliases=["Buchwald-Hartwig", "Pd-catalyzed amination"],
    category="C-N Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of amines with aryl halides to form aryl amines (C-N bonds).",
    conditions="Pd2(dba)3 or Pd(OAc)2 + bulky phosphine ligand (BINAP, XPhos, SPhos, RuPhos, BrettPhos, tBuXPhos), strong base (NaOtBu, Cs2CO3, K3PO4), toluene or dioxane, 80-110°C.",
    substrate="Primary/secondary amines, amides, sulfonamides + aryl bromides/chlorides/triflates. Ligand choice critical for selectivity.",
    limitations="Competitive β-hydride elimination with some substrates. Sensitive to steric effects. Ligand screening often necessary.",
    mechanism="Oxidative addition → amine coordination/deprotonation → reductive elimination.",
    refs="Guram, Buchwald, J. Am. Chem. Soc. 1994, 116, 7901. Louie, Hartwig, Tetrahedron Lett. 1995, 36, 3609.",
    related="Ullmann coupling, Chan-Lam coupling, Goldberg reaction",
    functional_groups="amine, aryl halide")

_nr("Hiyama Coupling",
    aliases=["Hiyama reaction"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of organosilanes with organic halides. Uses inexpensive, non-toxic silicon reagents.",
    conditions="Pd(0) catalyst, fluoride activator (TBAF, KF), THF or DMF, RT-80°C.",
    substrate="Aryl/vinylsilanes (trialkoxysilanes, silanols) + aryl halides.",
    limitations="Requires fluoride activation. Slower transmetalation than Sn, B, Zn. Competing protodesilylation.",
    mechanism="Oxidative addition → fluoride-assisted transmetalation → reductive elimination.",
    refs="Hatanaka, Hiyama, J. Org. Chem. 1988, 53, 918. Denmark, Sweis, Acc. Chem. Res. 2002, 35, 835.",
    related="Suzuki, Stille, Negishi",
    functional_groups="organosilane, aryl halide")

# --- C-C Bond Formation: Classical ---

_nr("Grignard Reaction",
    aliases=["Grignard addition"],
    category="C-C Bond Formation", type_="Addition",
    summary="Addition of organomagnesium halides (RMgX) to electrophiles (carbonyls, epoxides, CO2, etc.).",
    conditions="RMgX in THF or Et2O, add to electrophile at 0°C to -78°C, then warm. Requires anhydrous/inert atmosphere.",
    substrate="RMgX + aldehydes → 2° alcohols; + ketones → 3° alcohols; + esters → 3° alcohols (2 equiv); + CO2 → carboxylic acids; + epoxides → alcohols.",
    limitations="Incompatible with protic groups (OH, NH, COOH) and many electrophilic FGs. Cannot use protic solvents. Sensitive to moisture/air.",
    mechanism="Nucleophilic addition. Schlenk equilibrium: 2 RMgX ⇌ R2Mg + MgX2.",
    refs="Grignard, Compt. Rend. 1900, 130, 1322. Nobel Prize 1912.",
    related="Organolithium, Barbier, Reformatsky, Kumada coupling",
    functional_groups="organomagnesium, aldehyde, ketone, ester, epoxide")

_nr("Aldol Reaction",
    aliases=["Aldol condensation", "Aldol addition"],
    category="C-C Bond Formation", type_="Addition",
    summary="Addition of enolate or enol to aldehyde/ketone forming β-hydroxy carbonyl. Condensation gives α,β-unsaturated carbonyl.",
    conditions="Base-catalyzed: LDA, NaOH, KOH; acid-catalyzed: BF3·Et2O, TiCl4. Mukaiyama aldol: silyl enol ether + Lewis acid.",
    substrate="Enolizable carbonyl + aldehyde/ketone. Crossed aldol requires careful control (Evans auxiliary, Zimmerman-Traxler model).",
    limitations="Self-condensation. Mixtures of syn/anti diastereomers without chiral auxiliary. Retro-aldol under harsh conditions.",
    mechanism="Enolate formation → nucleophilic addition to carbonyl → β-hydroxy carbonyl. E1cb dehydration gives α,β-unsaturated product.",
    refs="Mukaiyama, Org. React. 1982, 28, 203. Evans, Aldrichimica Acta 1982, 15, 23.",
    related="Mukaiyama aldol, Claisen condensation, Reformatsky, Mannich",
    functional_groups="enolate, aldehyde, ketone")

_nr("Wittig Reaction",
    aliases=["Wittig olefination"],
    category="C-C Bond Formation", type_="Olefination",
    summary="Reaction of phosphonium ylides with aldehydes/ketones to form alkenes with loss of Ph3P=O.",
    conditions="Ylide from R3P+CH2R'X- + base (nBuLi, NaHMDS, KOtBu). Non-stabilized ylides: THF, -78°C. Stabilized ylides: RT.",
    substrate="Phosphonium salt + aldehyde (or reactive ketone). Selectivity: non-stabilized ylides → Z-alkene; stabilized ylides → E-alkene.",
    limitations="Removal of Ph3P=O can be difficult. E/Z selectivity not always predictable. Ketones react slowly. MW of phosphine oxide waste.",
    mechanism="Ylide + carbonyl → betaine → oxaphosphetane → [2+2] cycloreversion → alkene + Ph3P=O.",
    refs="Wittig, Geissler, Justus Liebigs Ann. Chem. 1953, 580, 44. Nobel Prize 1979. Maryanoff, Reitz, Chem. Rev. 1989, 89, 863.",
    related="HWE, Julia olefination, Peterson olefination, Tebbe olefination",
    functional_groups="phosphonium ylide, aldehyde, ketone")

_nr("Horner-Wadsworth-Emmons Reaction",
    aliases=["HWE", "HWE reaction", "Horner-Wadsworth-Emmons", "HWE olefination"],
    category="C-C Bond Formation", type_="Olefination",
    summary="Olefination using phosphonate-stabilized carbanions. More E-selective than Wittig; water-soluble byproduct.",
    conditions="Phosphonate ester + base (NaH, LiHMDS, DBU) in THF, 0°C to RT. Still-Gennari modification for Z-selectivity.",
    substrate="Phosphonate esters (Arbuzov product) + aldehydes. Excellent E-selectivity with standard conditions. Z-selective with Still-Gennari (CF3CH2O phosphonate + KHMDS/18-crown-6).",
    limitations="Less reactive with ketones. Requires preparation of phosphonate ester.",
    mechanism="Similar to Wittig but via open transition state. Phosphonate anion + aldehyde → β-hydroxyphosphonate → elimination.",
    refs="Horner, Hoffmann, Wippel, Chem. Ber. 1958, 91, 61. Wadsworth, Emmons, J. Am. Chem. Soc. 1961, 83, 1733. Still, Gennari, Tetrahedron Lett. 1983, 24, 4405.",
    related="Wittig, Julia olefination, Peterson",
    functional_groups="phosphonate ester, aldehyde")

_nr("Diels-Alder Reaction",
    aliases=["[4+2] cycloaddition", "Diels-Alder"],
    category="Cycloaddition", type_="Pericyclic",
    summary="[4+2] cycloaddition of conjugated diene with dienophile to form cyclohexene. Concerted, stereospecific.",
    conditions="Thermal. Lewis acid catalysis (AlCl3, BF3·Et2O, ZnCl2) accelerates and improves endo selectivity. Chiral Lewis acids for asymmetric DA.",
    substrate="Diene (s-cis conformation, electron-rich preferred) + dienophile (electron-poor preferred). Normal electron demand: EDG-diene + EWG-dienophile. Inverse ED possible.",
    limitations="Diene must adopt s-cis. Retro-DA at high temperature. Endo/exo selectivity. Hetero-DA possible with C=O, C=N dienophiles.",
    mechanism="Concerted [4πs + 2πs] pericyclic reaction. Suprafacial-suprafacial. Endo rule (Alder rule). Woodward-Hoffmann allowed thermally.",
    refs="Diels, Alder, Justus Liebigs Ann. Chem. 1928, 460, 98. Nobel Prize 1950. Nicolaou et al., Angew. Chem. Int. Ed. 2002, 41, 1668.",
    related="Hetero-Diels-Alder, 1,3-dipolar cycloaddition, [2+2] cycloaddition, retro-Diels-Alder",
    functional_groups="diene, dienophile")

_nr("Mannich Reaction",
    aliases=["Mannich condensation"],
    category="C-C Bond Formation", type_="Multicomponent",
    summary="Three-component reaction of aldehyde + amine + enolizable ketone to form β-amino carbonyl (Mannich base).",
    conditions="Acid catalysis (AcOH, HCl). Formaldehyde + secondary amine + ketone in water/ethanol. Organocatalytic asymmetric versions available.",
    substrate="Non-enolizable aldehyde (formaldehyde most common) + primary or secondary amine + enolizable carbonyl compound.",
    limitations="Mixtures with unsymmetrical ketones. Double Mannich possible. Retro-Mannich under some conditions.",
    mechanism="Iminium ion formation from aldehyde + amine → enol or enolate attacks iminium ion → β-amino carbonyl.",
    refs="Mannich, Krösche, Arch. Pharm. 1912, 250, 647. Arend et al., Angew. Chem. Int. Ed. 1998, 37, 1044.",
    related="Aldol, Strecker, Henry reaction",
    functional_groups="aldehyde, amine, enol")

_nr("Henry Reaction",
    aliases=["Nitroaldol reaction", "Henry nitroaldol"],
    category="C-C Bond Formation", type_="Addition",
    summary="Base-catalyzed addition of nitroalkane to aldehyde/ketone giving β-nitro alcohol.",
    conditions="Base (NaOH, K2CO3, DBU, guanidine). Asymmetric versions with Cu/BOX, rare earth catalysts.",
    substrate="Nitroalkane (nitromethane most common) + aldehyde. NO2 group can be transformed: → amine (reduction), → carbonyl (Nef), → alkene (elimination).",
    limitations="Retro-Henry possible. Elimination to nitroalkene can occur (useful or undesired). Mixtures of diastereomers.",
    mechanism="Deprotonation of nitroalkane → nitronate anion attacks carbonyl.",
    refs="Henry, Compt. Rend. 1895, 120, 1265. Luzzio, Tetrahedron 2001, 57, 915.",
    related="Aldol, Nitro-Mannich, Nef reaction",
    functional_groups="nitroalkane, aldehyde")

_nr("Reformatsky Reaction",
    aliases=["Reformatsky"],
    category="C-C Bond Formation", type_="Addition",
    summary="Zinc-mediated addition of α-haloesters to aldehydes/ketones giving β-hydroxy esters.",
    conditions="α-Bromoester + Zn dust in THF or benzene, reflux. Sonication or activated Zn (Rieke Zn) helpful.",
    substrate="α-Bromo ester + aldehyde/ketone. Better FG tolerance than Grignard due to lower reactivity of organozinc.",
    limitations="Less reactive than Grignard. Sometimes sluggish initiation. Diastereoselectivity moderate.",
    mechanism="Zn inserts into C-Br bond → organozinc enolate → addition to carbonyl.",
    refs="Reformatsky, Ber. Dtsch. Chem. Ges. 1887, 20, 1210. Rathke, Org. React. 1975, 22, 423.",
    related="Aldol, Grignard, Blaise reaction",
    functional_groups="α-haloester, aldehyde, ketone, zinc")

_nr("Baylis-Hillman Reaction",
    aliases=["Morita-Baylis-Hillman", "MBH reaction"],
    category="C-C Bond Formation", type_="Addition",
    summary="Coupling of activated alkene (α,β-unsaturated carbonyl) with aldehyde in presence of nucleophilic catalyst (DABCO, PPh3).",
    conditions="DABCO (5-100 mol%) or PPh3 in neat conditions or in solvent (DMSO, MeOH), RT. Often slow (days).",
    substrate="Acrylate, MVK, acrylonitrile + aldehydes. Product is α-methylene-β-hydroxy ester/ketone.",
    limitations="Very slow reaction rates (sometimes weeks). Can be accelerated with pressure or Lewis acid co-catalysts.",
    mechanism="Nucleophilic catalyst adds to alkene → zwitterionic enolate → aldol with aldehyde → proton transfer → elimination of catalyst.",
    refs="Baylis, Hillman, German Patent 2155113, 1972. Basavaiah et al., Chem. Rev. 2003, 103, 811.",
    related="Aldol, Rawal Diels-Alder",
    functional_groups="acrylate, aldehyde")

# --- Olefinations ---

_nr("Julia-Lythgoe Olefination",
    aliases=["Julia olefination", "Julia-Lythgoe", "Modified Julia", "Julia-Kocienski"],
    category="C-C Bond Formation", type_="Olefination",
    summary="Olefination via sulfone chemistry. Classic Julia: phenylsulfone + aldehyde → β-hydroxy sulfone → reductive elimination → E-alkene. Modified Julia-Kocienski: one-pot using benzothiazolyl (BT) or tetrazolyl (PT) sulfones.",
    conditions="Classic: PhSO2CHR + nBuLi + RCHO → β-hydroxy sulfone; then Na/Hg (or SmI2) + MeOH → E-alkene. Modified (Kocienski): BT/PT sulfone + KHMDS + RCHO, one pot, -78°C to RT.",
    substrate="Sulfone + aldehyde. E-selective (classic). Modified Julia: E or Z depending on sulfone type.",
    limitations="Classic requires two steps. Na(Hg) amalgam unpleasant. Modified Julia variable E/Z ratios.",
    mechanism="Classic: addition → β-elimination (reductive). Modified: addition → Smiles rearrangement → SO2 loss.",
    refs="Julia, Paris, Tetrahedron Lett. 1973, 14, 4833. Kocienski et al., Synlett 1998, 26.",
    related="Wittig, HWE, Peterson olefination",
    functional_groups="sulfone, aldehyde")

_nr("Peterson Olefination",
    aliases=["Peterson reaction"],
    category="C-C Bond Formation", type_="Olefination",
    summary="Olefination of aldehydes/ketones using α-silyl carbanions. Stereochemistry controlled by elimination conditions.",
    conditions="α-Silyl carbanion (from TMSCH2Li or TMSCH2MgCl) + carbonyl. Acidic workup → E-alkene. Basic workup → Z-alkene.",
    substrate="α-TMS organolithium/Grignard + aldehyde/ketone.",
    limitations="Requires preparation of α-silyl organometallic.",
    mechanism="Addition → β-hydroxy silane. Acid: anti-elimination (E). Base: syn-elimination via pentacoordinate Si (Z).",
    refs="Peterson, J. Org. Chem. 1968, 33, 780.",
    related="Wittig, HWE, Julia",
    functional_groups="α-silyl carbanion, aldehyde, ketone")

# --- Oxidations ---

_nr("Swern Oxidation",
    aliases=["Swern"],
    category="Oxidation", type_="Alcohol → Aldehyde/Ketone",
    summary="Oxidation of primary/secondary alcohols to aldehydes/ketones using DMSO/oxalyl chloride.",
    conditions="(COCl)2 + DMSO in CH2Cl2 at -78°C, then add alcohol, then Et3N. Warm to RT.",
    substrate="Primary → aldehyde (no over-oxidation). Secondary → ketone. Very mild, compatible with sensitive substrates.",
    limitations="Must maintain -78°C during DMSO activation (explosive side products above -60°C). Generates CO, CO2, Me2S (stinky). Stoichiometric reagents.",
    mechanism="DMSO activated by oxalyl chloride → chlorosulfonium salt. Alcohol attacks → alkoxysulfonium → Et3N-mediated intramolecular elimination → carbonyl + Me2S.",
    refs="Omura, Swern, Tetrahedron 1978, 34, 1651. Mancuso, Swern, Synthesis 1981, 165.",
    related="Dess-Martin, PCC, PDC, Ley (TPAP), Parikh-Doering, Pfitzner-Moffatt",
    functional_groups="alcohol")

_nr("Dess-Martin Periodinane",
    aliases=["Dess-Martin", "DMP oxidation", "Dess-Martin oxidation"],
    category="Oxidation", type_="Alcohol → Aldehyde/Ketone",
    summary="Mild oxidation of alcohols using Dess-Martin periodinane (DMP). Very popular for sensitive substrates.",
    conditions="1.1-1.5 equiv DMP in CH2Cl2, RT, 15-60 min. Add NaHCO3 for acid-sensitive substrates.",
    substrate="Primary → aldehyde. Secondary → ketone. Tolerates many FGs including alkenes, alkynes, silyl ethers.",
    limitations="Expensive reagent. Potentially explosive (contains I(V)). Shelf life limited. CH2Cl2 solvent usually required.",
    mechanism="Ligand exchange on hypervalent iodine → reductive elimination → carbonyl + reduced iodine species.",
    refs="Dess, Martin, J. Org. Chem. 1983, 48, 4155. Dess, Martin, J. Am. Chem. Soc. 1991, 113, 7277.",
    related="Swern, PCC, TPAP, IBX",
    functional_groups="alcohol")

_nr("Jones Oxidation",
    aliases=["Jones reagent", "CrO3/H2SO4"],
    category="Oxidation", type_="Alcohol → Carboxylic Acid/Ketone",
    summary="CrO3/H2SO4 in acetone. Oxidizes primary alcohols to carboxylic acids and secondary alcohols to ketones.",
    conditions="Jones reagent (CrO3 in dilute H2SO4) in acetone, 0°C to RT. Titrate to persistent orange color.",
    substrate="Primary alcohols → carboxylic acids. Secondary alcohols → ketones. Also oxidizes aldehydes to acids.",
    limitations="Strong conditions — over-oxidation of 1° alcohols to acids. Not compatible with acid-sensitive groups. Generates toxic Cr waste.",
    mechanism="Chromate ester formation → [2,3]-elimination.",
    refs="Bowden, Heilbron, Jones, J. Chem. Soc. 1946, 39.",
    related="PCC, PDC, Swern, Collins (CrO3·py2)",
    functional_groups="alcohol, aldehyde")

_nr("Sharpless Epoxidation",
    aliases=["Sharpless AE", "SAE", "asymmetric epoxidation"],
    category="Oxidation", type_="Asymmetric Epoxidation",
    summary="Ti-catalyzed asymmetric epoxidation of allylic alcohols using TBHP and diethyl tartrate (DET).",
    conditions="Ti(OiPr)4 (5-10 mol%), (+)- or (-)-DET or DIPT, TBHP, 4Å mol sieves, CH2Cl2, -20°C.",
    substrate="Allylic alcohols. Predictable face selectivity from mnemonic: draw allylic OH at bottom-right → (+)-DET attacks from top, (-)-DET from bottom.",
    limitations="Only works for allylic alcohols (OH directs). Non-allylic alkenes: use Jacobsen or Shi epoxidation. ee depends on alkene substitution pattern (trisubstituted best).",
    mechanism="Ti-tartrate complex forms chiral pocket → TBHP delivers oxygen to one face of alkene.",
    refs="Katsuki, Sharpless, J. Am. Chem. Soc. 1980, 102, 5974. Nobel Prize 2001.",
    related="Jacobsen epoxidation, Shi epoxidation, Sharpless dihydroxylation",
    functional_groups="allylic alcohol, alkene")

_nr("Sharpless Dihydroxylation",
    aliases=["Sharpless AD", "SAD", "asymmetric dihydroxylation", "AD-mix"],
    category="Oxidation", type_="Asymmetric Dihydroxylation",
    summary="Os-catalyzed asymmetric dihydroxylation of alkenes to syn-1,2-diols using AD-mix-α or AD-mix-β.",
    conditions="AD-mix-α (contains (DHQ)2PHAL + K2OsO4·2H2O + K3Fe(CN)6 + K2CO3) or AD-mix-β (contains (DHQD)2PHAL). tBuOH/H2O 1:1, 0°C.",
    substrate="Virtually any alkene. Mnemonic predicts face selectivity. Best ee with trans-disubstituted alkenes.",
    limitations="Os is toxic and expensive (catalytic amounts mitigate). Terminal alkenes give lower ee. MeSO2NH2 accelerates slow substrates.",
    mechanism="OsO4 + chiral amine ligand → [3+2] cycloaddition with alkene → osmate ester → hydrolysis → diol + Os(VI). Reoxidation by K3Fe(CN)6.",
    refs="Jacobsen, Marko, Mungall, Sharpless, J. Am. Chem. Soc. 1988, 110, 1968. Nobel Prize 2001.",
    related="Sharpless epoxidation, Upjohn dihydroxylation (OsO4/NMO), Sharpless aminohydroxylation",
    functional_groups="alkene")

_nr("Baeyer-Villiger Oxidation",
    aliases=["Baeyer-Villiger", "BV oxidation"],
    category="Oxidation", type_="Ketone → Ester/Lactone",
    summary="Oxidation of ketone to ester (or cyclic ketone to lactone) using peracid. Migratory aptitude determines regiochemistry.",
    conditions="mCPBA (1.1-3 equiv) in CH2Cl2, RT to reflux. Also: CF3CO3H, MMPP, H2O2/BF3.",
    substrate="Cyclic ketones → lactones. Acyclic ketones → esters. Migratory aptitude: 3° alkyl > cyclohexyl > 2° alkyl > phenyl > 1° alkyl > methyl.",
    limitations="Competing epoxidation if alkene present. mCPBA can be hazardous on scale. Acidic byproduct.",
    mechanism="Nucleophilic addition of peracid to ketone → Criegee intermediate → 1,2-migration with loss of carboxylate.",
    refs="Baeyer, Villiger, Ber. Dtsch. Chem. Ges. 1899, 32, 3625. Krow, Org. React. 1993, 43, 251.",
    related="mCPBA epoxidation, Beckmann rearrangement",
    functional_groups="ketone")

_nr("Pinnick Oxidation",
    aliases=["Pinnick", "Lindgren oxidation", "NaClO2 oxidation", "Kraus oxidation"],
    category="Oxidation", type_="Aldehyde → Carboxylic Acid",
    summary="Selective oxidation of aldehyde to carboxylic acid using NaClO2 (sodium chlorite). Does not affect alkenes, alkynes, or ketones.",
    conditions="NaClO2 (2-5 equiv), NaH2PO4 (buffer), 2-methyl-2-butene (HOCl scavenger), tBuOH/H2O, RT.",
    substrate="Aldehydes → carboxylic acids. Very chemoselective — alkenes, alcohols, ethers untouched.",
    limitations="2-Methyl-2-butene scavenger essential (otherwise HOCl causes side reactions). Chlorine gas can evolve.",
    mechanism="NaClO2 oxidizes aldehyde via chlorous acid. 2-Methyl-2-butene quenches HOCl byproduct.",
    refs="Bal, Pinnick, Tetrahedron 1981, 37, 2091. Lindgren, Nilsson, Acta Chem. Scand. 1973, 27, 888.",
    related="Jones, PDC, MnO2, TEMPO",
    functional_groups="aldehyde")

_nr("TEMPO Oxidation",
    aliases=["TEMPO", "Anelli oxidation", "TEMPO/bleach"],
    category="Oxidation", type_="Alcohol → Aldehyde",
    summary="Catalytic TEMPO with stoichiometric co-oxidant selectively oxidizes primary alcohols to aldehydes.",
    conditions="TEMPO (cat.), NaOCl (bleach) or PhI(OAc)2 or BAIB, CH2Cl2 or CH2Cl2/H2O biphasic. Zhao variant: TEMPO/BAIB/CH2Cl2.",
    substrate="Primary alcohols → aldehydes (no over-oxidation to acid). Secondary alcohols → ketones. Chemoselective for 1° over 2°.",
    limitations="Not effective for hindered secondary alcohols. Some substrates require different co-oxidant.",
    mechanism="Oxoammonium cation is active oxidant. TEMPO → oxoammonium → hydroxylamine (recycled by co-oxidant).",
    refs="Anelli et al., J. Org. Chem. 1987, 52, 2559. De Mico et al., J. Org. Chem. 1997, 62, 6974.",
    related="Swern, Dess-Martin, PCC, Ley (TPAP)",
    functional_groups="primary alcohol")

_nr("Ley Oxidation",
    aliases=["TPAP", "TPAP/NMO", "Ley-Griffith"],
    category="Oxidation", type_="Alcohol → Aldehyde/Ketone",
    summary="Catalytic TPAP (Pr4N+ RuO4-) with NMO as co-oxidant. Mild, selective, good for sensitive substrates.",
    conditions="TPAP (2-5 mol%), NMO (1.5 equiv), 4Å mol sieves, CH2Cl2, RT.",
    substrate="Primary → aldehyde. Secondary → ketone. Excellent chemoselectivity.",
    limitations="Sensitive to scale (exothermic). Molecular sieves essential. Can be sluggish with hindered substrates.",
    mechanism="Ru(VII) → Ru(V) during oxidation, regenerated by NMO.",
    refs="Ley et al., Synthesis 1994, 639. Griffith et al., J. Chem. Soc. Chem. Commun. 1987, 1625.",
    related="Swern, Dess-Martin, TEMPO",
    functional_groups="alcohol")

_nr("Ozonolysis",
    aliases=["Ozonolysis", "Criegee mechanism"],
    category="Oxidation", type_="Alkene Cleavage",
    summary="Cleavage of C=C double bonds with ozone to give aldehydes/ketones (reductive workup) or carboxylic acids (oxidative workup).",
    conditions="O3 bubbled through substrate in CH2Cl2 or MeOH at -78°C until blue color persists. Reductive workup: Me2S or PPh3. Oxidative workup: H2O2.",
    substrate="Any alkene. Product depends on substitution pattern and workup. Reductive: aldehydes/ketones. Oxidative: carboxylic acids/ketones.",
    limitations="O3 is toxic and generated on-site. Must quench peroxidic intermediates. Scale-up safety concerns.",
    mechanism="[3+2] cycloaddition → molozonide → retro-[3+2] → carbonyl + carbonyl oxide → [3+2] → ozonide. Workup cleaves ozonide.",
    refs="Criegee, Angew. Chem. Int. Ed. 1975, 14, 745.",
    related="Lemieux-Johnson (OsO4/NaIO4), dihydroxylation/periodate",
    functional_groups="alkene")

_nr("Wacker Oxidation",
    aliases=["Wacker process", "Wacker-Tsuji"],
    category="Oxidation", type_="Alkene → Ketone",
    summary="Pd(II)-catalyzed oxidation of terminal alkenes to methyl ketones (Markovnikov).",
    conditions="PdCl2 (10 mol%), CuCl (co-oxidant), O2, DMF/H2O, RT-60°C.",
    substrate="Terminal alkenes → methyl ketones. Anti-Markovnikov selectivity possible with specific ligands.",
    limitations="Internal alkenes unreactive. Can give aldehyde byproduct. Over-oxidation possible.",
    mechanism="Pd(II) coordinates alkene → nucleophilic attack of water (Markovnikov) → β-hydride elimination → enol → ketone. Pd(0) reoxidized by CuCl2/O2.",
    refs="Smidt et al., Angew. Chem. 1959, 71, 176. Tsuji, Synthesis 1984, 369.",
    related="Heck, hydroboration-oxidation",
    functional_groups="terminal alkene")

# --- Reductions ---

_nr("Birch Reduction",
    aliases=["Birch"],
    category="Reduction", type_="Aromatic → 1,4-Cyclohexadiene",
    summary="Dissolving metal reduction of aromatic rings to 1,4-cyclohexadienes using Na/Li in NH3(l) + alcohol.",
    conditions="Na or Li metal in liquid NH3 (-33°C) with tBuOH or iPrOH as proton source. THF co-solvent.",
    substrate="Arenes. EDG (OMe, NR2): 2,5-diene (unconjugated). EWG (COOH, COOR): 1,4-diene (conjugated with EWG).",
    limitations="Requires liquid NH3 (hazardous). Over-reduction possible. Na/NH3 incompatible with many FGs (halides, epoxides).",
    mechanism="Solvated electrons reduce aromatic ring → radical anion → protonation → radical → second electron → carbanion → protonation.",
    refs="Birch, J. Chem. Soc. 1944, 430. Rabideau, Marcinow, Org. React. 1992, 42, 1.",
    related="Benkeser reduction, catalytic hydrogenation",
    functional_groups="arene")

_nr("Wolff-Kishner Reduction",
    aliases=["Wolff-Kishner", "Huang-Minlon modification"],
    category="Reduction", type_="Carbonyl → Methylene",
    summary="Reduction of aldehydes/ketones to alkanes via hydrazone formation and base-mediated N2 loss.",
    conditions="Hydrazine (NH2NH2) + KOH in ethylene glycol (or diethylene glycol), reflux (200°C). Huang-Minlon: one-pot.",
    substrate="Aldehydes/ketones → CH2. Complementary to Clemmensen (acidic conditions).",
    limitations="High temperatures. Base-sensitive groups not tolerated. Long reaction times. Not for acid-sensitive substrates (use this instead of Clemmensen).",
    mechanism="Hydrazone formation → deprotonation by base → carbanion → loss of N2 → protonation → alkane.",
    refs="Wolff, Justus Liebigs Ann. Chem. 1912, 394, 86. Kishner, J. Russ. Phys. Chem. Soc. 1911, 43, 582. Huang-Minlon, J. Am. Chem. Soc. 1946, 68, 2487.",
    related="Clemmensen reduction, thioacetal/Raney Ni desulfurization",
    functional_groups="aldehyde, ketone")

_nr("Clemmensen Reduction",
    aliases=["Clemmensen"],
    category="Reduction", type_="Carbonyl → Methylene",
    summary="Reduction of aldehydes/ketones to alkanes using Zn(Hg)/conc. HCl. Acidic complement to Wolff-Kishner.",
    conditions="Zn(Hg) amalgam, concentrated HCl, reflux.",
    substrate="Aryl ketones → aryl alkanes. Particularly useful for Friedel-Crafts acylation products.",
    limitations="Strong acid conditions. Not suitable for acid-sensitive substrates. Base-sensitive substrates OK (complement to W-K).",
    mechanism="Heterogeneous — exact mechanism debated. Occurs on zinc surface. Not via alcohol intermediate.",
    refs="Clemmensen, Ber. Dtsch. Chem. Ges. 1913, 46, 1837.",
    related="Wolff-Kishner, thioacetal reduction",
    functional_groups="aldehyde, ketone")

_nr("Luche Reduction",
    aliases=["Luche", "NaBH4/CeCl3"],
    category="Reduction", type_="Selective 1,2-Reduction",
    summary="CeCl3-mediated NaBH4 reduction. Selectively reduces ketones (1,2-addition) in presence of conjugated enones.",
    conditions="CeCl3·7H2O (1.1 equiv) + NaBH4 (0.5-1 equiv) in MeOH, 0°C to RT.",
    substrate="α,β-Unsaturated ketones → allylic alcohols (1,2-reduction). Highly selective over 1,4-reduction.",
    limitations="MeOH solvent required. CeCl3 must be dry for some applications.",
    mechanism="Ce(III) activates carbonyl toward 1,2-attack. NaBH4 provides hydride. Hard-hard interaction favored.",
    refs="Luche, J. Am. Chem. Soc. 1978, 100, 2226. Gemal, Luche, J. Am. Chem. Soc. 1981, 103, 5454.",
    related="NaBH4, L-Selectride, CBS reduction, DIBAL",
    functional_groups="enone, ketone")

_nr("CBS Reduction",
    aliases=["Corey-Bakshi-Shibata", "CBS", "oxazaborolidine"],
    category="Reduction", type_="Asymmetric Ketone Reduction",
    summary="Enantioselective reduction of prochiral ketones using chiral oxazaborolidine catalyst and BH3.",
    conditions="(R) or (S)-CBS catalyst (5-20 mol%) + BH3·THF or BH3·SMe2, toluene or THF, -20°C to RT.",
    substrate="Prochiral ketones → chiral secondary alcohols. Predictable stereochemistry via transition state model.",
    limitations="Requires stoichiometric BH3. Sensitive to moisture. Catalyst can be expensive.",
    mechanism="Oxazaborolidine activates both BH3 (Lewis acid on B) and ketone → six-membered TS → enantioselective hydride delivery.",
    refs="Corey, Bakshi, Shibata, J. Am. Chem. Soc. 1987, 109, 5551.",
    related="Luche, Noyori hydrogenation, Alpine-Borane, DIP-Chloride",
    functional_groups="ketone")

# --- Rearrangements ---

_nr("Claisen Rearrangement",
    aliases=["Claisen", "oxy-Cope", "Ireland-Claisen", "Johnson-Claisen", "Eschenmoser-Claisen"],
    category="Rearrangement", type_="[3,3]-Sigmatropic",
    summary="[3,3]-Sigmatropic rearrangement of allyl vinyl ethers to γ,δ-unsaturated carbonyls. Many variants.",
    conditions="Thermal: 150-250°C (aromatic Claisen) or RT-80°C (aliphatic with activation). Ireland-Claisen: LDA/TMSCl then warm. Johnson-Claisen: CH3C(OEt)3, cat. acid.",
    substrate="Allyl vinyl ether → γ,δ-unsaturated carbonyl. Chair-like TS gives predictable stereochemistry (E-enolate → anti product).",
    limitations="High temperatures for unactivated substrates. Requires allyl vinyl ether synthesis.",
    mechanism="Concerted [3,3]-sigmatropic shift through chair-like transition state. Suprafacial.",
    refs="Claisen, Ber. Dtsch. Chem. Ges. 1912, 45, 3157. Ireland, Mueller, Willard, J. Am. Chem. Soc. 1976, 98, 2868.",
    related="Cope, oxy-Cope, [2,3]-Wittig, Overman",
    functional_groups="allyl vinyl ether")

_nr("Beckmann Rearrangement",
    aliases=["Beckmann"],
    category="Rearrangement", type_="Oxime → Amide",
    summary="Rearrangement of oximes to amides (or lactams from cyclic ketoximes) under acidic conditions.",
    conditions="Acid catalyst: H2SO4, PCl5, SOCl2, TsCl/base, Beckmann mixture (AcOH/HCl/Ac2O). Heat.",
    substrate="Ketoximes → amides. Cyclic ketoximes → lactams (e.g., cyclohexanone oxime → caprolactam for nylon-6). Anti-periplanar group migrates.",
    limitations="Regioselectivity: anti group migrates (E/Z of oxime determines product). Strong acid conditions.",
    mechanism="Protonation of OH → loss of water → migration of anti group → nitrilium ion → water attack → amide.",
    refs="Beckmann, Ber. Dtsch. Chem. Ges. 1886, 19, 988. Gawley, Org. React. 1988, 35, 1.",
    related="Baeyer-Villiger, Curtius, Hofmann, Schmidt",
    functional_groups="oxime, ketone")

_nr("Curtius Rearrangement",
    aliases=["Curtius degradation"],
    category="Rearrangement", type_="Acyl Azide → Isocyanate",
    summary="Thermal rearrangement of acyl azides to isocyanates (which can be trapped as amines, carbamates, or ureas).",
    conditions="Acyl azide heated in inert solvent (toluene, 80-110°C). Modern: DPPA + Et3N generates acyl azide in situ. Trap with ROH → carbamate, H2O → amine, RNH2 → urea.",
    substrate="Carboxylic acids → amines (one carbon shorter). Via: RCOOH → RCOCl → RCON3 → RNCO → RNH2.",
    limitations="Acyl azides are potentially explosive. DPPA method much safer. HN3 generation possible.",
    mechanism="Concerted 1,2-shift with loss of N2 from acyl azide → isocyanate.",
    refs="Curtius, Ber. Dtsch. Chem. Ges. 1890, 23, 3023. Shioiri, Ninomiya, Yamada, J. Am. Chem. Soc. 1972, 94, 6203 (DPPA).",
    related="Hofmann, Lossen, Schmidt, Wolff rearrangement",
    functional_groups="acyl azide, carboxylic acid")

# --- Key Named Reactions ---

_nr("Mitsunobu Reaction",
    aliases=["Mitsunobu"],
    category="Substitution", type_="SN2 with Inversion",
    summary="Converts ROH + HNu → RNu with inversion (SN2) using PPh3/DIAD (or DEAD). Converts alcohols to esters, ethers, amines, azides, etc.",
    conditions="PPh3 (1.1 equiv), DIAD or DEAD (1.1 equiv), nucleophile (pKa < 11 for pronucleophile), THF, 0°C to RT.",
    substrate="Secondary alcohols → inverted product. Primary alcohols work but no stereo consequence. Nucleophiles: carboxylic acids, phenols, phthalimide, HN3, sulfonamides.",
    limitations="pKa requirement for pronucleophile (< 11-13). Stoichiometric PPh3/DIAD waste (difficult to remove). Triphenylphosphine oxide byproduct.",
    mechanism="PPh3 + DIAD → betaine → activates alcohol (oxyphosphonium) → SN2 by nucleophile → inversion.",
    refs="Mitsunobu, Yamada, Bull. Chem. Soc. Jpn. 1967, 40, 2380. Mitsunobu, Synthesis 1981, 1.",
    related="Appel reaction, Williamson ether synthesis, Gabriel synthesis",
    functional_groups="alcohol, carboxylic acid, phenol")

_nr("Appel Reaction",
    aliases=["Appel"],
    category="Substitution", type_="OH → Halide",
    summary="Conversion of alcohols to alkyl halides using PPh3/CX4 (X = Cl, Br). Inversion of stereochemistry.",
    conditions="PPh3 + CBr4 (→ RBr) or CCl4 (→ RCl) in CH2Cl2, 0°C to RT.",
    substrate="Primary and secondary alcohols → alkyl halides with inversion.",
    limitations="Stoichiometric PPh3 waste. Ph3P=O removal. Side reactions with tertiary alcohols.",
    mechanism="PPh3 + CX4 → [Ph3PX]+ + CHX3- → alcohol attacks P → oxyphosphonium → SN2 by halide.",
    refs="Appel, Angew. Chem. Int. Ed. 1975, 14, 801.",
    related="Mitsunobu, HBr, SOCl2, PBr3",
    functional_groups="alcohol")

_nr("Williamson Ether Synthesis",
    aliases=["Williamson"],
    category="Substitution", type_="SN2 Ether Formation",
    summary="SN2 reaction of alkoxide with alkyl halide to form ethers.",
    conditions="ROH + NaH (or KOH) → RO⁻; then add R'X (primary halide/tosylate), THF or DMF, 0°C to RT.",
    substrate="Alkoxide + primary (or methyl) alkyl halide. Secondary halides → elimination competes.",
    limitations="SN2 mechanism: primary halides work best. Secondary/tertiary → E2 elimination. Requires strong base to form alkoxide.",
    mechanism="Deprotonation → alkoxide → SN2 on alkyl halide → ether.",
    refs="Williamson, J. Chem. Soc. 1852, 4, 229.",
    related="Mitsunobu, Ullmann ether synthesis",
    functional_groups="alcohol, alkyl halide")

_nr("Gabriel Synthesis",
    aliases=["Gabriel"],
    category="Substitution", type_="Primary Amine Synthesis",
    summary="Synthesis of primary amines from alkyl halides using potassium phthalimide, avoiding over-alkylation.",
    conditions="Potassium phthalimide + RX (primary halide) in DMF, RT-80°C → N-alkylphthalimide → hydrazinolysis (N2H4, Ing-Manske) or acid hydrolysis → primary amine.",
    substrate="Primary alkyl halides → primary amines. Avoids polyalkylation problem of direct NH3 alkylation.",
    limitations="Only for primary amines. SN2 conditions (no tertiary/neopentyl halides). Hydrazinolysis step required.",
    mechanism="SN2 of phthalimide anion on alkyl halide. Deprotection by hydrazinolysis or acid.",
    refs="Gabriel, Ber. Dtsch. Chem. Ges. 1887, 20, 2224.",
    related="Staudinger reduction (azide → amine), reductive amination, Curtius",
    functional_groups="alkyl halide, amine")

_nr("Olefin Metathesis",
    aliases=["Grubbs metathesis", "RCM", "ring-closing metathesis", "cross metathesis", "CM", "ROMP"],
    category="C-C Bond Formation", type_="Metathesis",
    summary="Redistribution of alkene fragments catalyzed by Ru or Mo carbene complexes. Includes RCM, CM, and ROMP.",
    conditions="Grubbs I (PCy3)2Cl2Ru=CHPh, Grubbs II (IMes)(PCy3)Cl2Ru=CHPh, Hoveyda-Grubbs II. CH2Cl2 or toluene, RT-40°C, often under N2.",
    substrate="RCM: dienes → cyclic alkenes (5-7 membered rings best). CM: two alkenes → cross product. ROMP: cyclic olefins → polymers.",
    limitations="E/Z selectivity in CM. Ru contamination of products. Catalyst poisoning by amines, thiols, phosphines. Ethylene removal drives equilibrium.",
    mechanism="Chauvin mechanism: [2+2] cycloaddition/cycloreversion between metal carbene and alkene → metallacyclobutane → products.",
    refs="Grubbs, Handbook of Metathesis, 2003. Nobel Prize 2005 (Grubbs, Schrock, Chauvin).",
    related="Wittig, HWE, Julia",
    functional_groups="alkene")

_nr("Click Chemistry",
    aliases=["CuAAC", "Huisgen cycloaddition", "azide-alkyne cycloaddition", "click reaction"],
    category="Cycloaddition", type_="1,3-Dipolar Cycloaddition",
    summary="Cu(I)-catalyzed azide-alkyne cycloaddition (CuAAC) to form 1,4-disubstituted 1,2,3-triazoles regioselectively.",
    conditions="CuSO4·5H2O (5-10 mol%) + sodium ascorbate (reductant), tBuOH/H2O, RT. Or: CuI + DIPEA.",
    substrate="Organic azide + terminal alkyne → 1,4-triazole. SPAAC (strain-promoted): no Cu needed with cyclooctynes.",
    limitations="Requires terminal alkyne (internal alkynes unreactive). Cu toxicity in biological applications (use SPAAC). Azide safety (small organic azides can be explosive).",
    mechanism="Cu(I) forms copper acetylide → coordinates azide → stepwise cycloaddition → regioselective 1,4-triazole.",
    refs="Rostovtsev, Green, Fokin, Sharpless, Angew. Chem. Int. Ed. 2002, 41, 2596. Meldal, Tornøe, Chem. Rev. 2008, 108, 2952.",
    related="Huisgen 1,3-dipolar cycloaddition, Diels-Alder",
    functional_groups="azide, terminal alkyne, triazole")

_nr("Strecker Synthesis",
    aliases=["Strecker"],
    category="C-C Bond Formation", type_="Multicomponent",
    summary="Three-component synthesis of α-amino acids from aldehyde + amine (or NH3) + HCN. Gives α-aminonitrile → hydrolysis → amino acid.",
    conditions="RCHO + NH4Cl + NaCN in H2O, RT. Or: RCHO + amine + TMSCN + Lewis acid catalyst.",
    substrate="Aldehydes + amines + cyanide source → α-aminonitriles → α-amino acids (after hydrolysis).",
    limitations="HCN toxicity. Racemic (asymmetric Strecker with chiral catalysts available). Hydrolysis step needed.",
    mechanism="Aldehyde + amine → imine → cyanide addition to imine → α-aminonitrile.",
    refs="Strecker, Justus Liebigs Ann. Chem. 1850, 75, 27.",
    related="Mannich, Ugi, Passerini",
    functional_groups="aldehyde, amine, cyanide")

_nr("Ugi Reaction",
    aliases=["Ugi four-component reaction", "Ugi 4CR"],
    category="C-C Bond Formation", type_="Multicomponent",
    summary="Four-component condensation of aldehyde + amine + carboxylic acid + isocyanide → α-acylamino amide.",
    conditions="Mix all four components in MeOH, RT, 24-48h. High atom economy.",
    substrate="Aldehyde + primary amine + carboxylic acid + isocyanide. Enormous diversity — combinatorial chemistry workhorse.",
    limitations="Diastereoselectivity control difficult. Product often requires further cyclization/functionalization.",
    mechanism="Imine formation → protonation by acid → α-addition of isocyanide → acyl migration (Mumm rearrangement) → product.",
    refs="Ugi, Steinbrückner, Angew. Chem. 1960, 72, 267. Dömling, Ugi, Angew. Chem. Int. Ed. 2000, 39, 3168.",
    related="Strecker, Passerini, Mannich, Biginelli",
    functional_groups="aldehyde, amine, carboxylic acid, isocyanide")

_nr("Sandmeyer Reaction",
    aliases=["Sandmeyer"],
    category="Substitution", type_="Diazonium → Halide/CN/etc.",
    summary="Conversion of aryl diazonium salts to aryl halides (Cl, Br, I), cyanides, or other groups via Cu(I) catalysis.",
    conditions="ArNH2 + NaNO2/HCl (0°C) → ArN2+. Then: CuCl → ArCl; CuBr → ArBr; CuCN → ArCN; KI → ArI (no Cu needed).",
    substrate="Aryl amines → aryl halides, nitriles, and other derivatives via diazonium chemistry.",
    limitations="Diazonium salts can be explosive (especially dry). Must keep cold. Yields sometimes moderate.",
    mechanism="Single electron transfer from Cu(I) → aryl radical + N2 + Cu(II) → radical rebound with halide.",
    refs="Sandmeyer, Ber. Dtsch. Chem. Ges. 1884, 17, 1633.",
    related="Balz-Schiemann (→ ArF), Meerwein arylation, Schiemann",
    functional_groups="aryl amine, diazonium salt")

_nr("Fischer Esterification",
    aliases=["Fischer ester synthesis", "acid-catalyzed esterification"],
    category="Functional Group Interconversion", type_="Acid + Alcohol → Ester",
    summary="Acid-catalyzed esterification of carboxylic acid with alcohol. Equilibrium reaction — drive with excess alcohol or remove water.",
    conditions="RCOOH + R'OH, cat. H2SO4 or p-TsOH, reflux with Dean-Stark trap or molecular sieves.",
    substrate="Carboxylic acids + primary/secondary alcohols → esters. Equilibrium favored by excess of one reactant or water removal.",
    limitations="Equilibrium — reversible. Slow with bulky substrates. Tertiary alcohols → elimination instead.",
    mechanism="Protonation of C=O → nucleophilic addition of alcohol → proton transfer → loss of water → ester.",
    refs="Fischer, Speier, Ber. Dtsch. Chem. Ges. 1895, 28, 3252.",
    related="Steglich, Yamaguchi, Mukaiyama, DCC coupling",
    functional_groups="carboxylic acid, alcohol")

_nr("Steglich Esterification",
    aliases=["DCC coupling", "Steglich"],
    category="Functional Group Interconversion", type_="Acid + Alcohol → Ester",
    summary="DCC (or EDC)-mediated esterification with DMAP catalyst. Mild conditions for sensitive substrates.",
    conditions="RCOOH + R'OH + DCC (1.1 equiv) + DMAP (cat., 5-20 mol%), CH2Cl2, 0°C to RT.",
    substrate="Carboxylic acids + alcohols (including hindered). Also for amide bond formation (peptide coupling with HOBt).",
    limitations="Dicyclohexylurea (DCU) byproduct can be hard to remove. Racemization possible with amino acids (add HOBt). EDC gives water-soluble urea (easier removal).",
    mechanism="DCC activates acid → O-acylisourea → DMAP catalysis via acylpyridinium → nucleophilic substitution by alcohol.",
    refs="Neises, Steglich, Angew. Chem. Int. Ed. 1978, 17, 522.",
    related="Fischer, Yamaguchi, Mukaiyama esterification, amide coupling",
    functional_groups="carboxylic acid, alcohol")

_nr("Corey-Fuchs Reaction",
    aliases=["Corey-Fuchs"],
    category="C-C Bond Formation", type_="Aldehyde → Alkyne",
    summary="Two-step conversion of aldehyde to terminal alkyne via gem-dibromo olefin intermediate.",
    conditions="Step 1: RCHO + CBr4 + PPh3 → RCH=CBr2. Step 2: 2 equiv nBuLi, THF, -78°C → RC≡CH (or RC≡CLi for trapping).",
    substrate="Aldehydes → terminal alkynes (one carbon homologation).",
    limitations="Stoichiometric PPh3/CBr4 in step 1. Strong base (nBuLi) in step 2.",
    mechanism="Step 1: Wittig-type. Step 2: deprotonation + α-elimination → vinylidene carbene → rearrangement → alkyne.",
    refs="Corey, Fuchs, Tetrahedron Lett. 1972, 13, 3769.",
    related="Seyferth-Gilbert/Ohira-Bestmann, Wittig",
    functional_groups="aldehyde, alkyne")

_nr("Ohira-Bestmann Reagent",
    aliases=["Ohira-Bestmann", "Seyferth-Gilbert homologation", "dimethyl-1-diazo-2-oxopropylphosphonate"],
    category="C-C Bond Formation", type_="Aldehyde → Alkyne",
    summary="One-pot conversion of aldehyde to terminal alkyne using dimethyl-1-diazo-2-oxopropylphosphonate + K2CO3/MeOH.",
    conditions="Ohira-Bestmann reagent + K2CO3, MeOH, RT, 12-24h. Much milder than Corey-Fuchs.",
    substrate="Aldehydes → terminal alkynes. Tolerates many FGs.",
    limitations="Reagent preparation (from dimethyl-2-oxopropylphosphonate + TsN3). Only works with aldehydes (not ketones).",
    mechanism="HWE/Wittig-like → diazo intermediate → Wolff rearrangement → vinylidene carbene → 1,2-H shift → alkyne.",
    refs="Ohira, Synth. Commun. 1989, 19, 561. Müller, Liepold, Bestmann, Synlett 1996, 521.",
    related="Corey-Fuchs, Seyferth-Gilbert",
    functional_groups="aldehyde, alkyne")

_nr("Staudinger Reduction",
    aliases=["Staudinger"],
    category="Reduction", type_="Azide → Amine",
    summary="Reduction of organic azides to primary amines using triphenylphosphine.",
    conditions="RN3 + PPh3, THF/H2O, RT. Iminophosphorane intermediate hydrolyzed to amine.",
    substrate="Alkyl/aryl azides → primary amines. Staudinger ligation (no hydrolysis): used in bioconjugation.",
    limitations="Stoichiometric PPh3 (Ph3P=O waste). Azide must be accessible.",
    mechanism="PPh3 attacks terminal N of azide → phosphazide → loss of N2 → iminophosphorane → hydrolysis → amine + Ph3P=O.",
    refs="Staudinger, Meyer, Helv. Chim. Acta 1919, 2, 635.",
    related="Hydrogenation (Pd/C, H2), Gabriel synthesis, Curtius",
    functional_groups="azide, amine")

_nr("Arndt-Eistert Homologation",
    aliases=["Arndt-Eistert"],
    category="Homologation", type_="Acid → Homologated Acid",
    summary="One-carbon homologation of carboxylic acids via diazoketone and Wolff rearrangement.",
    conditions="RCOOH → RCOCl (SOCl2) → RCOCHN2 (CH2N2) → Ag2O/H2O or hv → RCH2COOH.",
    substrate="Carboxylic acids → one-carbon-homologated acids (or esters, amides if trapped with ROH, RNH2).",
    limitations="Diazomethane is toxic, explosive, carcinogenic. Silver catalysis needed. TMS-diazomethane as safer alternative.",
    mechanism="Wolff rearrangement: diazoketone → ketocarbene → ketene (1,2-shift) → trapping by nucleophile.",
    refs="Arndt, Eistert, Ber. Dtsch. Chem. Ges. 1935, 68, 200.",
    related="Corey-Fuchs, Curtius, Wolff rearrangement",
    functional_groups="carboxylic acid, diazoketone")

_nr("Finkelstein Reaction",
    aliases=["Finkelstein"],
    category="Substitution", type_="Halide Exchange",
    summary="SN2 halide exchange: RCl/RBr + NaI → RI + NaCl/NaBr. Driven by precipitation of NaCl/NaBr in acetone.",
    conditions="RX + NaI (excess) in acetone, reflux. NaCl/NaBr precipitate (insoluble in acetone) drives equilibrium.",
    substrate="Primary alkyl chlorides/bromides → iodides. Classic example of Le Chatelier's principle.",
    limitations="Only works well in acetone (NaCl/NaBr insoluble). Primarily for primary substrates (SN2).",
    mechanism="SN2: iodide attacks carbon, halide leaves.",
    refs="Finkelstein, Ber. Dtsch. Chem. Ges. 1910, 43, 1528.",
    related="Appel, Williamson",
    functional_groups="alkyl halide")

_nr("Reductive Amination",
    aliases=["reductive amination"],
    category="C-N Bond Formation", type_="Carbonyl + Amine → Amine",
    summary="Formation of C-N bond by condensation of carbonyl with amine followed by in situ reduction of imine/iminium.",
    conditions="RCHO or R2CO + R'NH2 → imine → NaBH3CN or NaBH(OAc)3 (pH 6-7) → amine. One-pot or stepwise.",
    substrate="Aldehydes (more reactive) or ketones + primary/secondary amines → secondary/tertiary amines.",
    limitations="NaBH3CN: toxic (HCN generation in acid). NaBH(OAc)3: milder, preferred. Competing carbonyl reduction. Ketones slower than aldehydes.",
    mechanism="Condensation → imine (or iminium ion) → hydride reduction → amine.",
    refs="Borch, Bernstein, Durst, J. Am. Chem. Soc. 1971, 93, 2897. Abdel-Magid et al., J. Org. Chem. 1996, 61, 3849.",
    related="Gabriel, Leuckart-Wallach, Eschweiler-Clarke, Mannich",
    functional_groups="aldehyde, ketone, amine")

_nr("Friedel-Crafts Alkylation",
    aliases=["Friedel-Crafts alkylation", "FC alkylation"],
    category="C-C Bond Formation", type_="Electrophilic Aromatic Substitution",
    summary="Lewis acid-catalyzed alkylation of aromatic rings with alkyl halides, alkenes, or alcohols.",
    conditions="ArH + RCl + AlCl3 (>1 equiv) in CH2Cl2 or CS2, 0°C to RT. Other Lewis acids: FeCl3, BF3·Et2O, ZnCl2.",
    substrate="Electron-rich arenes + alkyl halides/alkenes. Product is activated → polyalkylation common.",
    limitations="Over-alkylation (product more reactive than starting material). Carbocation rearrangements. Deactivated arenes unreactive. Not for arenes with NH2, NHR, OH (they coordinate to Lewis acid).",
    mechanism="Lewis acid generates carbocation from RX → electrophilic aromatic substitution → arenium ion → deprotonation.",
    refs="Friedel, Crafts, Compt. Rend. 1877, 84, 1392.",
    related="Friedel-Crafts acylation, Blanc chloromethylation",
    functional_groups="arene, alkyl halide")

_nr("Friedel-Crafts Acylation",
    aliases=["Friedel-Crafts acylation", "FC acylation"],
    category="C-C Bond Formation", type_="Electrophilic Aromatic Substitution",
    summary="Lewis acid-mediated acylation of aromatic rings with acyl halides or anhydrides. No rearrangement or over-acylation.",
    conditions="ArH + RCOCl + AlCl3 (>1 equiv), CH2Cl2, 0°C to RT. AlCl3 is consumed (stoichiometric — complexes with ketone product).",
    substrate="Electron-rich arenes + acyl chlorides/anhydrides → aryl ketones.",
    limitations="Requires >1 equiv AlCl3 (product-Lewis acid complex). Deactivated arenes unreactive. Works for ketone synthesis where alkylation fails (no rearrangement, no polyacylation).",
    mechanism="AlCl3 + RCOCl → acylium ion (or complex) → EAS → arenium ion → deprotonation.",
    refs="Friedel, Crafts, J. Chem. Soc. 1877, 32, 725.",
    related="Friedel-Crafts alkylation, Fries rearrangement, Haworth synthesis",
    functional_groups="arene, acyl halide")

_nr("Wittig Rearrangement",
    aliases=["[1,2]-Wittig", "[2,3]-Wittig"],
    category="Rearrangement", type_="[1,2] or [2,3]-Sigmatropic",
    summary="Base-induced rearrangement of ethers. [1,2]-Wittig: α-alkoxy carbanion → 1,2-shift → alkoxide. [2,3]-Wittig: α-allyloxy carbanion → [2,3]-sigmatropic → homoallylic alkoxide.",
    conditions="nBuLi or LDA, THF, -78°C. Substrate must have acidic proton α to oxygen.",
    substrate="[1,2]: R-O-CHR' (with stabilizing group) → alcohol. [2,3]: allyl ethers with α-proton → homoallylic alcohols.",
    limitations="Requires α-proton on ether. [1,2] is stepwise (radical pair). [2,3] is concerted and stereospecific.",
    mechanism="[1,2]: homolytic C-O cleavage → radical pair → recombination. [2,3]: concerted suprafacial [2,3]-sigmatropic shift.",
    refs="Wittig, Löhmann, Justus Liebigs Ann. Chem. 1942, 550, 260. Nakai, Mikami, Chem. Rev. 1986, 86, 885.",
    related="Claisen rearrangement, Stevens rearrangement",
    functional_groups="ether, allyl ether")

_nr("Swern Oxidation",
    aliases=["Swern"],
    category="Oxidation", type_="Alcohol → Aldehyde/Ketone",
    summary="Oxidation of primary/secondary alcohols to aldehydes/ketones using DMSO/oxalyl chloride.",
    conditions="(COCl)2 + DMSO in CH2Cl2 at -78°C, then add alcohol, then Et3N. Warm to RT.",
    substrate="Primary → aldehyde. Secondary → ketone.",
    limitations="Must maintain -78°C. Generates Me2S (smell).",
    mechanism="DMSO activated by oxalyl chloride → alkoxysulfonium → elimination.",
    refs="Omura, Swern, Tetrahedron 1978, 34, 1651.",
    related="Dess-Martin, PCC, TEMPO, Parikh-Doering",
    functional_groups="alcohol")

_nr("Oppenauer Oxidation",
    aliases=["Oppenauer"],
    category="Oxidation", type_="Alcohol → Ketone (Al-mediated)",
    summary="Aluminum alkoxide-catalyzed oxidation of secondary alcohols to ketones using acetone as hydride acceptor. Reverse of MPV reduction.",
    conditions="Al(OiPr)3, excess acetone, reflux.",
    substrate="Secondary alcohols → ketones. Mild — tolerates alkenes, alkynes.",
    limitations="Equilibrium reaction. Excess acetone needed. Not practical for primary alcohols.",
    mechanism="Meerwein-Ponndorf-Verley in reverse: Al coordinates both substrate and acetone → hydride transfer → ketone + isopropanol.",
    refs="Oppenauer, Recl. Trav. Chim. Pays-Bas 1937, 56, 137.",
    related="MPV reduction, Ley (TPAP), Swern",
    functional_groups="secondary alcohol")

_nr("Meerwein-Ponndorf-Verley Reduction",
    aliases=["MPV", "Meerwein-Ponndorf-Verley", "MPV reduction"],
    category="Reduction", type_="Ketone → Alcohol (Al-mediated)",
    summary="Aluminum alkoxide-catalyzed transfer hydrogenation of ketones using isopropanol as hydride source.",
    conditions="Al(OiPr)3, excess iPrOH, reflux.",
    substrate="Ketones → secondary alcohols. Highly chemoselective for C=O over C=C.",
    limitations="Equilibrium. Does not reduce esters, amides, etc.",
    mechanism="Six-membered cyclic TS: Al coordinates ketone substrate and iPrOH → hydride transfer → alcohol + acetone.",
    refs="Meerwein, Schmidt, Justus Liebigs Ann. Chem. 1925, 444, 221.",
    related="Oppenauer oxidation, Noyori transfer hydrogenation",
    functional_groups="ketone")

_nr("Barton Decarboxylation",
    aliases=["Barton", "Barton radical decarboxylation"],
    category="Radical", type_="Decarboxylation",
    summary="Radical decarboxylation of carboxylic acids via Barton esters (thiohydroxamate esters). Very general radical generation method.",
    conditions="RCOOH → mixed anhydride with Barton's reagent (N-hydroxy-2-thiopyridone sodium salt) → photolysis or AIBN → radical R· + CO2. Trap with Bu3SnH or other radical trap.",
    substrate="Carboxylic acids → R· radical → R-H, R-Br, R-SPy, etc. Works with 1°, 2°, 3° acids.",
    limitations="Requires Barton ester preparation. Radical conditions (sensitivity to O2). Tin reagents toxic.",
    mechanism="Barton ester → photolysis → carboxyl radical → β-scission (loss of CO2) → carbon radical → trapping.",
    refs="Barton, Crich, Motherwell, J. Chem. Soc. Chem. Commun. 1983, 939.",
    related="Hunsdiecker, Kolbe electrolysis",
    functional_groups="carboxylic acid")

_nr("Fischer Indole Synthesis",
    aliases=["Fischer indole"],
    category="Heterocycle Formation", type_="Indole Synthesis",
    summary="Synthesis of indoles from aryl hydrazones via [3,3]-sigmatropic rearrangement and loss of NH3.",
    conditions="Aryl hydrazone + Lewis or Brønsted acid (ZnCl2, HCl, AcOH, polyphosphoric acid), heat (100-200°C).",
    substrate="Aryl hydrazine + aldehyde/ketone → hydrazone → indole. Regiochemistry predictable for unsymmetrical ketones.",
    limitations="High temperatures. Regiochemistry issues with unsymmetrical ketones (migrate more substituted carbon). NH3 loss.",
    mechanism="Tautomerization to ene-hydrazine → [3,3]-sigmatropic rearrangement → re-aromatization → loss of NH3 → indole.",
    refs="Fischer, Jourdan, Ber. Dtsch. Chem. Ges. 1883, 16, 2241.",
    related="Bartoli, Leimgruber-Batcho, Madelung indole synthesis",
    functional_groups="aryl hydrazine, aldehyde, ketone, indole")

_nr("Paal-Knorr Synthesis",
    aliases=["Paal-Knorr furan", "Paal-Knorr pyrrole", "Paal-Knorr"],
    category="Heterocycle Formation", type_="Furan/Pyrrole/Thiophene Synthesis",
    summary="Synthesis of furans, pyrroles, or thiophenes from 1,4-dicarbonyls. Choice of reagent determines product.",
    conditions="1,4-Diketone + acid catalyst → furan. + Primary amine → pyrrole. + Lawesson's reagent or P4S10 → thiophene.",
    substrate="1,4-Dicarbonyl compounds → 5-membered heterocycles.",
    limitations="Requires accessible 1,4-dicarbonyl (can be synthesized via Stetter, 1,4-addition, etc.).",
    mechanism="Double condensation/cyclodehydration. Furan: acid-catalyzed enolization → ring closure → dehydration × 2.",
    refs="Paal, Ber. Dtsch. Chem. Ges. 1884, 17, 2756. Knorr, Ber. Dtsch. Chem. Ges. 1884, 17, 1635.",
    related="Hantzsch pyrrole, Feist-Bénary furan",
    functional_groups="1,4-diketone, furan, pyrrole, thiophene")

# --- Additional Cross-Coupling ---

_nr("Buchwald-Hartwig Amination (C-O)",
    aliases=["Buchwald etherification", "C-O coupling"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed C-O bond formation between aryl halides and alcohols/phenols.",
    conditions="Pd2(dba)3 or Pd(OAc)2, bulky phosphine (RuPhos, tBuXPhos), Cs2CO3 or NaOtBu, toluene, 80-110°C.",
    substrate="Aryl halides + phenols, primary/secondary alcohols.",
    limitations="Less developed than C-N variant. Competing β-hydride elimination with 2° alcohols.",
    mechanism="Oxidative addition → coordination/deprotonation of alcohol → reductive elimination.",
    refs="Vorogushin, Huang, Buchwald, JACS 2005, 127, 8146.",
    related="Buchwald-Hartwig, Ullmann, Chan-Lam",
    functional_groups="aryl halide, alcohol, ether")

_nr("Chan-Lam Coupling",
    aliases=["Chan-Lam-Evans", "copper-mediated coupling"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Cu-mediated coupling of boronic acids with amines, phenols, or thiols under mild conditions.",
    conditions="Cu(OAc)2 (stoich. or cat.), Et3N or pyridine, CH2Cl2 or DMF, RT, air atmosphere. No Pd required.",
    substrate="ArB(OH)2 + amines, phenols, amides, or alcohols. Open-flask conditions (O2 as terminal oxidant).",
    limitations="Often requires stoichiometric Cu. Homo-coupling side reaction. Less reliable than Pd-catalyzed methods.",
    mechanism="Transmetalation → Cu(II)/Cu(III) oxidation → reductive elimination. Cu(I)/Cu(II)/Cu(III) cycle.",
    refs="Chan, Lam, JACS 1998, 120, 8657. Evans, Katz, West, Tetrahedron Lett. 1998, 39, 2937.",
    related="Buchwald-Hartwig, Ullmann, Suzuki",
    functional_groups="boronic acid, amine, phenol")

_nr("Ullmann Coupling",
    aliases=["Ullmann reaction", "Ullmann-Goldberg"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Cu-catalyzed coupling of aryl halides with amines, phenols, or thiols. Modern variants use mild conditions.",
    conditions="CuI (5-20 mol%), diamine ligand (DMEDA, proline), K2CO3 or Cs2CO3, DMSO or dioxane, 80-110°C.",
    substrate="Aryl iodides/bromides + amines, phenols, amides, thiols. Original: 2 ArX → Ar-Ar (high T, Cu powder).",
    limitations="Less reliable than Pd catalysis. Requires iodides typically. High temperatures.",
    mechanism="Cu(I)/Cu(III) oxidative addition/reductive elimination cycle (debated: SET vs concerted).",
    refs="Ullmann, Ber. 1903, 36, 2382. Ma, Cai, Acc. Chem. Res. 2008, 41, 1450.",
    related="Buchwald-Hartwig, Chan-Lam",
    functional_groups="aryl halide, amine, phenol")

_nr("Tsuji-Trost Reaction",
    aliases=["Trost allylation", "Pd-allyl", "allylic substitution"],
    category="C-C Bond Formation", type_="Substitution",
    summary="Pd-catalyzed allylic substitution. Nucleophilic attack on π-allyl Pd complex.",
    conditions="Pd(PPh3)4 or [Pd(allyl)Cl]2, ligand (PPh3, DPPF, or chiral PHOX for asymmetric), THF or CH2Cl2, RT-60°C.",
    substrate="Allylic acetates, carbonates, halides + soft nucleophiles (malonates, amines) or hard nucleophiles (organozinc).",
    limitations="Regioselectivity depends on nucleophile hardness (soft → less substituted, hard → more substituted). Linear/branched selectivity.",
    mechanism="Oxidative addition → π-allyl complex formation → nucleophilic attack (outer sphere for soft, inner sphere for hard nuc).",
    refs="Trost, Fullerton, JACS 1973, 95, 292. Trost, Van Vranken, Chem. Rev. 1996, 96, 395.",
    related="Heck, Suzuki",
    functional_groups="allylic acetate, malonate, amine")

_nr("Hiyama Coupling",
    aliases=["Hiyama cross-coupling", "Hiyama-Denmark"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of organosilanes with organic halides. Requires fluoride activation.",
    conditions="Pd(0) cat., TBAF or CsF (activator), THF or DMF, 50-100°C.",
    substrate="RSiR'3 (trialkoxysilanes, silanols) + ArX. Organosilanes are cheap, stable, low toxicity.",
    limitations="Requires activator (F⁻ or base) for transmetalation. Lower reactivity than boronic acids.",
    mechanism="Oxidative addition → F⁻-activated transmetalation (pentacoordinate Si) → reductive elimination.",
    refs="Hiyama, Hatanaka, Pure Appl. Chem. 1994, 66, 1471. Denmark, Regens, Acc. Chem. Res. 2008, 41, 1486.",
    related="Suzuki, Stille, Negishi",
    functional_groups="organosilane, aryl halide")

_nr("Liebeskind-Srogl Coupling",
    aliases=["thioester coupling"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of thioester C-S bond with boronic acids under neutral conditions (no base needed).",
    conditions="Pd(PPh3)4, CuTC (Cu(I) thienylcarboxylate), THF, RT-50°C. pH-neutral.",
    substrate="Thioester + ArB(OH)2. Produces ketone. Also works with thioethers.",
    limitations="Requires stoichiometric Cu co-catalyst. Thioester substrate preparation.",
    mechanism="Oxidative addition of C-S bond → Cu-assisted transmetalation → reductive elimination.",
    refs="Liebeskind, Srogl, JACS 2000, 122, 11260.",
    related="Suzuki, Stille, Fukuyama coupling",
    functional_groups="thioester, boronic acid, ketone")

_nr("Fukuyama Coupling",
    aliases=["Fukuyama reduction"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of thioesters with organozinc reagents to form ketones. Selective mono-addition.",
    conditions="Pd(PPh3)4, organozinc (RZnX), THF, 0°C to RT.",
    substrate="Thioesters + RZnX → ketones. Does NOT over-add to give alcohol (unlike aldehydes + RZnX).",
    limitations="Requires preparation of organozinc reagent. Thioester substrate required.",
    mechanism="Oxidative addition of C-S bond → transmetalation with RZnX → reductive elimination.",
    refs="Tokuyama, Yokoshima, Mori, Li, Fukuyama, JACS 1998, 120, 11235.",
    related="Negishi, Weinreb amide, Liebeskind-Srogl",
    functional_groups="thioester, organozinc, ketone")

# --- Olefination Reactions ---

_nr("Still-Gennari Olefination",
    aliases=["Still-Gennari"],
    category="C-C Bond Formation", type_="Olefination",
    summary="Modified HWE reaction that gives Z-alkenes selectively using bis(trifluoroethyl) phosphonate.",
    conditions="(CF3CH2O)2P(O)CH2CO2R, KHMDS, 18-crown-6, THF, -78°C.",
    substrate="Aldehydes → Z-α,β-unsaturated esters (>95:5 Z:E).",
    limitations="Expensive reagent. Requires low temperature and crown ether. Limited to ester-stabilized ylides.",
    mechanism="Modified HWE. Kinetic control at low T. cis-Oxaphosphetane forms preferentially due to CF3 groups.",
    refs="Still, Gennari, Tetrahedron Lett. 1983, 24, 4405.",
    related="HWE, Wittig, Julia, Ando",
    functional_groups="aldehyde, phosphonate, Z-alkene")

_nr("Takai Olefination",
    aliases=["Takai reaction", "Takai-Utimoto"],
    category="C-C Bond Formation", type_="Olefination",
    summary="CrCl2-mediated conversion of aldehydes to E-vinyl halides using CHX3.",
    conditions="CrCl2 (6 equiv), CHI3 or CHBr3, THF, 0°C to RT. cat. NiCl2 sometimes added.",
    substrate="RCHO → RCH=CHX (X = I or Br). Highly E-selective.",
    limitations="Stoichiometric Cr (toxic, but Cr(II) less so). Over-reduction possible.",
    mechanism="Cr(II) reduces CHX3 to Cr-carbenoid → Wittig-like addition to aldehyde.",
    refs="Takai, Nitta, Utimoto, JACS 1986, 108, 7408.",
    related="Wittig, HWE, Nozaki-Hiyama-Kishi",
    functional_groups="aldehyde, vinyl halide")

_nr("Tebbe Olefination",
    aliases=["Tebbe reagent", "Petasis olefination"],
    category="C-C Bond Formation", type_="Olefination",
    summary="Ti-mediated methylenation of carbonyls (including esters, amides) to terminal alkenes using Tebbe or Petasis reagent.",
    conditions="Tebbe reagent (Cp2TiCl2/AlMe3 complex) or Petasis reagent (Cp2TiMe2), THF, -40°C to RT, pyridine base.",
    substrate="Ketones, esters, lactones, amides → methylene (C=CH2). Works on carbonyls that Wittig cannot.",
    limitations="Only methylenation (no other alkylidenes with Tebbe). Sensitive reagent. Stoichiometric Ti.",
    mechanism="α-Elimination → Cp2Ti=CH2 (titanium carbene) → [2+2] cycloaddition with C=O → retro-[2+2] fragmentation.",
    refs="Tebbe, Parshall, Reddy, JACS 1978, 100, 3611. Petasis, Bzowej, JACS 1990, 112, 6392.",
    related="Wittig, Lombardo, Nysted",
    functional_groups="ketone, ester, alkene")

# --- Oxidation Reactions ---

_nr("Corey-Kim Oxidation",
    aliases=["Corey-Kim"],
    category="Oxidation", type_="Alcohol Oxidation",
    summary="Oxidation of alcohols using NCS/DMS, followed by Et3N. Mild alternative to Swern.",
    conditions="NCS, DMS, CH2Cl2, -25°C → add alcohol → Et3N, warm to RT.",
    substrate="1° alcohols → aldehydes. 2° alcohols → ketones.",
    limitations="Similar to Swern. DMS byproduct (odor). NCS must be fresh.",
    mechanism="NCS activates DMS → sulfonium salt → alkoxysulfonium ylide → intramolecular elimination.",
    refs="Corey, Kim, JACS 1972, 94, 7586.",
    related="Swern, Pfitzner-Moffatt, Parikh-Doering",
    functional_groups="alcohol, aldehyde, ketone")

_nr("Parikh-Doering Oxidation",
    aliases=["Doering oxidation", "SO3-pyridine oxidation"],
    category="Oxidation", type_="Alcohol Oxidation",
    summary="Mild DMSO-based oxidation using SO3·pyridine complex. Room temperature, no cryogenic conditions needed.",
    conditions="SO3·py, DMSO, Et3N, CH2Cl2, 0°C to RT.",
    substrate="1° → aldehyde, 2° → ketone. Very mild conditions.",
    limitations="DMSO required (difficult to remove). Slower than Swern.",
    mechanism="SO3·py activates DMSO → alkoxysulfonium intermediate → Kornblum elimination via Et3N.",
    refs="Parikh, Doering, JACS 1967, 89, 5505.",
    related="Swern, Corey-Kim, Pfitzner-Moffatt",
    functional_groups="alcohol, aldehyde, ketone")

_nr("IBX Oxidation",
    aliases=["IBX", "2-iodoxybenzoic acid oxidation"],
    category="Oxidation", type_="Alcohol Oxidation",
    summary="Mild, selective alcohol oxidation using IBX (2-iodoxybenzoic acid). Works in DMSO.",
    conditions="IBX (1.1-2 equiv), DMSO, RT, 1-12h. Filtered to remove iodosobenzoic acid byproduct.",
    substrate="1° → aldehyde (no over-oxidation!), 2° → ketone. Tolerates many functional groups.",
    limitations="IBX insoluble in most solvents except DMSO. Explosive when heated dry. Use freshly prepared.",
    mechanism="Hypervalent iodine(V). Ligand exchange at I(V) → [2,3]-sigmatropic rearrangement → elimination.",
    refs="Frigerio, Santagostino, J. Org. Chem. 1999, 64, 4537.",
    related="Dess-Martin, Swern, TEMPO",
    functional_groups="alcohol, aldehyde, ketone")

_nr("Rubottom Oxidation",
    aliases=["Rubottom"],
    category="Oxidation", type_="α-Hydroxylation",
    summary="α-Hydroxylation of enol silyl ethers using mCPBA to give α-hydroxy ketones.",
    conditions="Silyl enol ether + mCPBA, CH2Cl2, 0°C to RT → α-silyloxy ketone. Then TBAF or acid for deprotection.",
    substrate="Silyl enol ethers → α-hydroxy ketones. Regiospecific (defined by enolization).",
    limitations="Requires preparation of silyl enol ether (regiochemistry must be set). Over-oxidation possible.",
    mechanism="Electrophilic epoxidation of silyl enol ether → epoxide → 1,4-O→O silyl migration → α-silyloxyketone.",
    refs="Rubottom, Vazquez, Pelegrina, Tetrahedron Lett. 1974, 15, 4319.",
    related="Davis oxaziridine, Sharpless AD",
    functional_groups="silyl enol ether, α-hydroxy ketone")

_nr("Lemieux-Johnson Oxidation",
    aliases=["Lemieux-Johnson", "OsO4-NaIO4 cleavage"],
    category="Oxidation", type_="Alkene Cleavage",
    summary="One-step oxidative cleavage of alkenes to carbonyl compounds using OsO4/NaIO4.",
    conditions="cat. OsO4 (1-5 mol%), NaIO4 (2-3 equiv), dioxane/H2O or THF/H2O, RT.",
    substrate="Alkenes → 2 carbonyl fragments (aldehydes and/or ketones). Cleaves C=C completely.",
    limitations="Stoichiometric NaIO4 needed. Over-oxidation of aldehydes possible. OsO4 is toxic and expensive.",
    mechanism="OsO4 dihydroxylation → diol → NaIO4 cleaves diol → 2 carbonyls + regenerates OsO4.",
    refs="Lemieux, Johnson, J. Org. Chem. 1956, 21, 478.",
    related="Ozonolysis, Sharpless dihydroxylation",
    functional_groups="alkene, aldehyde, ketone")

_nr("Shi Epoxidation",
    aliases=["Shi asymmetric epoxidation"],
    category="Oxidation", type_="Asymmetric Epoxidation",
    summary="Organocatalytic asymmetric epoxidation of trans-disubstituted and trisubstituted alkenes using fructose-derived ketone.",
    conditions="Shi catalyst (30 mol%), Oxone, K2CO3, Na2B4O7 buffer, CH3CN/H2O, 0°C.",
    substrate="trans-Disubstituted, trisubstituted alkenes. ee >90% typically.",
    limitations="Less effective for cis- and terminal alkenes. High catalyst loading. pH-sensitive.",
    mechanism="Ketone + Oxone → dioxirane (chiral) → concerted O-transfer to alkene.",
    refs="Tu, Wang, Frohn, He, Yu, Tang, Shi, JACS 2003, 125, 14354.",
    related="Sharpless epoxidation, Jacobsen epoxidation, mCPBA",
    functional_groups="alkene, epoxide")

_nr("Jacobsen Epoxidation",
    aliases=["Jacobsen-Katsuki"],
    category="Oxidation", type_="Asymmetric Epoxidation",
    summary="Mn(III)-salen-catalyzed asymmetric epoxidation of unfunctionalized cis-disubstituted alkenes.",
    conditions="Mn(salen) cat. (1-5 mol%), NaOCl or mCPBA (oxidant), CH2Cl2, 0°C to RT.",
    substrate="cis-Disubstituted and trisubstituted alkenes (complementary to Sharpless AE which requires allylic alcohol).",
    limitations="trans-Alkenes give lower ee. Terminal alkenes variable. Catalyst decomposition over time.",
    mechanism="Mn(III) → Mn(V)=O (via oxidant) → concerted (though stepwise debated) O-transfer to alkene. Salen controls facial selectivity.",
    refs="Zhang, Loebach, Wilson, Jacobsen, JACS 1990, 112, 2801.",
    related="Sharpless epoxidation, Shi epoxidation, mCPBA",
    functional_groups="alkene, epoxide")

_nr("Upjohn Dihydroxylation",
    aliases=["OsO4 dihydroxylation", "osmylation"],
    category="Oxidation", type_="Dihydroxylation",
    summary="Catalytic OsO4-mediated syn-dihydroxylation of alkenes using NMO as co-oxidant.",
    conditions="cat. OsO4 (1-5 mol%), NMO (1.5 equiv), acetone/H2O or tBuOH/H2O, RT.",
    substrate="Any alkene → syn-1,2-diol. Electron-rich alkenes react faster. Stereospecific syn addition.",
    limitations="OsO4 is toxic and volatile. Can over-oxidize to cleave diol if NaIO4 is present.",
    mechanism="[3+2] cycloaddition of OsO4 across alkene → osmate(VI) ester → hydrolysis → diol + Os(VI). NMO reoxidizes Os(VI) to OsO4.",
    refs="Van Rheenen, Kelly, Cha, Tetrahedron Lett. 1976, 17, 1973.",
    related="Sharpless AD, Lemieux-Johnson",
    functional_groups="alkene, diol")

_nr("Wohl-Ziegler Bromination",
    aliases=["NBS bromination", "allylic bromination"],
    category="Oxidation", type_="Halogenation",
    summary="Allylic or benzylic bromination using NBS (N-bromosuccinimide) via radical mechanism.",
    conditions="NBS (1 equiv), radical initiator (AIBN, BPO, or light), CCl4 or CH2Cl2, reflux or hν.",
    substrate="Allylic, benzylic C-H → C-Br. Maintains alkene geometry. Regioselective for weakest C-H.",
    limitations="Mixture of products if multiple allylic positions. Competing addition to alkene if [Br2] too high. Over-bromination.",
    mechanism="Radical chain: initiation → Br• abstracts allylic/benzylic H → C• + NBS → C-Br + succinimidyl radical → propagation.",
    refs="Wohl, Ber. 1919, 52, 51. Ziegler, Angew. Chem. 1942, 55, 339.",
    related="radical bromination, Appel",
    functional_groups="allylic C-H, benzylic C-H, bromide")

# --- Reduction Reactions ---

_nr("Staudinger Ligation",
    aliases=["Staudinger reaction (amine)"],
    category="Reduction", type_="Azide Reduction",
    summary="Reduction of organic azides to amines using PPh3, via an iminophosphorane intermediate.",
    conditions="PPh3, THF/H2O, RT. Hydrolyze the iminophosphorane with water.",
    substrate="RN3 → RNH2. Very selective, tolerates many functional groups.",
    limitations="Produces Ph3P=O (stoichiometric, difficult to remove). The Bertozzi variant traps the intermediate for bioconjugation.",
    mechanism="PPh3 attacks azide → phosphazide → loss of N2 → iminophosphorane (R-N=PPh3) → hydrolysis → RNH2 + O=PPh3.",
    refs="Staudinger, Meyer, Helv. Chim. Acta 1919, 2, 635.",
    related="Curtius, Gabriel synthesis, hydrogenation",
    functional_groups="azide, amine, phosphine")

_nr("Noyori Asymmetric Hydrogenation",
    aliases=["BINAP hydrogenation", "Noyori hydrogenation"],
    category="Reduction", type_="Asymmetric Hydrogenation",
    summary="Ru-BINAP-catalyzed asymmetric hydrogenation of functionalized alkenes and ketones.",
    conditions="Ru(OAc)2(BINAP) or RuCl2(BINAP), H2 (50-100 atm), MeOH or iPrOH, 25-50°C.",
    substrate="β-Keto esters → β-hydroxy esters. Functionalized alkenes (enamides, allylic alcohols). ee >95% typical.",
    limitations="Requires coordinating group near C=C or C=O for high ee. High H2 pressure. Expensive catalyst.",
    mechanism="Substrate coordinates via functional group → H2 oxidative addition → migratory insertion (enantioselective face).",
    refs="Noyori, Ohkuma, Kitamura, JACS 1987, 109, 5856. Nobel Prize 2001.",
    related="CBS reduction, Corey-Bakshi-Shibata, Rh-DuPhos",
    functional_groups="ketone, alkene, β-hydroxy ester")

_nr("Bouveault-Blanc Reduction",
    aliases=["Bouveault-Blanc", "Na/alcohol reduction"],
    category="Reduction", type_="Ester Reduction",
    summary="Classical reduction of esters to alcohols using Na metal in an alcohol solvent (replaced by LiAlH4 in modern chemistry).",
    conditions="Na metal, EtOH or nBuOH, reflux. Historical method.",
    substrate="Esters → primary alcohols + alcohol from ester alkyl group.",
    limitations="Dangerous (Na + protic solvent). Superseded by LiAlH4 and DIBAL-H.",
    mechanism="Single-electron transfer from Na → radical anion → further reduction → alkoxide → protonation.",
    refs="Bouveault, Blanc, Compt. Rend. 1903, 136, 1676.",
    related="LiAlH4, DIBAL-H, Birch reduction",
    functional_groups="ester, alcohol")

_nr("Midland Reduction",
    aliases=["Alpine-Borane", "Midland alpine borane"],
    category="Reduction", type_="Asymmetric Reduction",
    summary="Asymmetric reduction of prochiral ketones using B-3-pinanyl-9-BBN (Alpine-Borane).",
    conditions="Alpine-Borane (from 9-BBN + α-pinene), THF, 0°C to RT. Slow, but high ee for propargylic ketones.",
    substrate="α,β-Ynones, α-halo ketones → chiral propargylic/secondary alcohols. Best for substrates with small/large differentiation.",
    limitations="Slow reaction. Limited substrate scope (best for acetylenic ketones). Requires large excess sometimes.",
    mechanism="Intramolecular hydride delivery via 6-membered TS. Boat-like transition state determines enantioselectivity.",
    refs="Midland, Greer, Tramontano, Zderic, JACS 1979, 101, 2352.",
    related="CBS reduction, Noyori, DIP-Cl",
    functional_groups="ketone, propargylic alcohol")

_nr("Meerwein Reduction",
    aliases=["Meerwein-Ponndorf-Verley", "MPV reduction"],
    category="Reduction", type_="Carbonyl Reduction",
    summary="Al(OiPr)3-catalyzed transfer hydrogenation of ketones/aldehydes using iPrOH as hydride source.",
    conditions="Al(OiPr)3, iPrOH, reflux. Also: Ln(OiPr)3, Sm(OiPr)3, or Ru catalysts for modern variants.",
    substrate="Ketones → 2° alcohols, aldehydes → 1° alcohols. Chemoselective: does not reduce esters, alkenes, etc.",
    limitations="Equilibrium reaction (remove acetone to drive forward). Slow with Al(OiPr)3.",
    mechanism="6-membered cyclic TS: Al coordinates both substrate C=O and iPrO-H → hydride transfer via Zimmerman-Traxler-like TS.",
    refs="Meerwein, Schmidt, Annalen 1925, 444, 221.",
    related="Oppenauer (reverse), NaBH4, CBS",
    functional_groups="ketone, aldehyde, alcohol")

# ---- v7.0 additions: 65+ new reactions ----

_nr("Bamford-Stevens Reaction",
    aliases=["Bamford-Stevens", "Shapiro Reaction"],
    category="C-C Bond Formation", type_="Elimination",
    summary="Tosylhydrazones → alkenes via diazo intermediates under basic conditions.",
    conditions="Tosylhydrazone + NaOMe or NaH, heat (150–200 °C) or nBuLi (Shapiro variant, milder). Shapiro: 2 eq nBuLi, –78 °C → rt.",
    substrate="Ketone-derived tosylhydrazones → alkenes. Shapiro gives less-substituted alkene (kinetic).",
    limitations="Bamford-Stevens: carbene rearrangements at high T. Shapiro requires strong base.",
    mechanism="Base removes NH → diazo intermediate → loss of N₂ → carbene (protic) or vinyl anion (Shapiro: base abstracts second H).",
    refs="Bamford, Stevens, J. Chem. Soc. 1952, 4735. Shapiro, Heath, JACS 1967, 89, 5734.",
    related="Wolff-Kishner, Julia Olefination",
    functional_groups="ketone, alkene, tosylhydrazone")

_nr("Bischler-Napieralski Reaction",
    aliases=["Bischler-Napieralski"],
    category="Heterocycle Formation", type_="Ring Closure",
    summary="Cyclodehydration of β-arylethylamides to 3,4-dihydroisoquinolines.",
    conditions="POCl3 or P2O5, reflux in toluene/xylene. Modern: Tf2O, 2-fluoropyridine.",
    substrate="β-Arylethylamides → dihydroisoquinolines. Electron-rich arenes cyclize faster.",
    limitations="Electron-poor arenes fail. Requires acyl group on nitrogen.",
    mechanism="POCl3 activates carbonyl → electrophilic cyclization → dehydration → imine.",
    refs="Bischler, Napieralski, Ber. 1893, 26, 1903.",
    related="Pictet-Spengler, Pomeranz-Fritsch",
    functional_groups="amide, isoquinoline")

_nr("Buchner Ring Expansion",
    aliases=["Buchner reaction"],
    category="Rearrangement", type_="Ring Expansion",
    summary="Carbene insertion into benzene ring → cycloheptatriene (tropylium precursor).",
    conditions="α-Diazo esters + Rh2(OAc)4, or Cu catalysis. Also photolytic.",
    substrate="Benzene, substituted arenes → cycloheptatrienes/norcaradienes.",
    limitations="Mixtures of isomers. Competing C-H insertion.",
    mechanism="Rh-carbenoid cyclopropanates benzene → norcaradiene ↔ cycloheptatriene electrocyclic equilibrium.",
    refs="Buchner, Curtius, Ber. 1885, 18, 2377.",
    related="Diels-Alder, Claisen Rearrangement",
    functional_groups="arene, carbene, cycloheptatriene")

_nr("Cannizzaro Reaction",
    aliases=["Cannizzaro"],
    category="Oxidation", type_="Disproportionation",
    summary="Base-mediated disproportionation of non-enolizable aldehydes → alcohol + carboxylate.",
    conditions="Conc. NaOH or KOH (50%), rt or mild heat. Crossed Cannizzaro: HCHO as reductant.",
    substrate="Aldehydes without α-H (ArCHO, R3CCHO, HCHO). Crossed variant uses formaldehyde.",
    limitations="Only non-enolizable aldehydes. 50% max yield per product (disproportionation).",
    mechanism="OH⁻ adds to C=O → tetrahedral alkoxide → hydride transfer to second aldehyde → alcohol + carboxylate.",
    refs="Cannizzaro, Liebigs Ann. 1853, 88, 129.",
    related="Tischenko, Aldol",
    functional_groups="aldehyde, carboxylic acid, alcohol")

_nr("Chichibabin Amination",
    aliases=["Chichibabin reaction"],
    category="C-N Bond Formation", type_="Nucleophilic Aromatic Substitution",
    summary="Direct amination of pyridines and related heterocycles with NaNH2.",
    conditions="NaNH2, toluene or mineral oil, 100–110 °C. Also: LiNH2, KNH2.",
    substrate="Pyridine → 2-aminopyridine. Also quinoline, isoquinoline, pyrimidine.",
    limitations="Strong base required. May give mixtures with polysubstituted substrates.",
    mechanism="NaNH2 attacks C-2 of pyridine → addition → NaH eliminated → 2-aminopyridine.",
    refs="Chichibabin, Zeide, JRCS 1914, 46, 1216.",
    related="Buchwald-Hartwig, Gabriel",
    functional_groups="pyridine, amine")

_nr("Claisen Condensation",
    aliases=["Claisen condensation"],
    category="C-C Bond Formation", type_="Condensation",
    summary="Base-mediated self-condensation of esters to give β-keto esters.",
    conditions="NaOEt/EtOH or LDA/THF, –78 °C → rt. Crossed: use LDA for controlled selectivity.",
    substrate="Esters with α-H → β-keto esters. Dieckmann: intramolecular (5- or 6-membered rings).",
    limitations="Requires α-H. Self-condensation gives mixtures with crossed variants.",
    mechanism="Base deprotonates α-H → enolate attacks ester C=O → tetrahedral intermediate → alkoxide leaves → β-keto ester.",
    refs="Claisen, Ber. 1887, 20, 655.",
    related="Aldol, Knoevenagel, Dieckmann",
    functional_groups="ester, β-keto ester")

_nr("Comins' Reagent",
    aliases=["Comins reagent", "vinyl triflate formation"],
    category="Functional Group Interconversion", type_="Triflation",
    summary="N-(5-Chloro-2-pyridyl)triflimide converts enolates to vinyl triflates.",
    conditions="LDA, THF, –78 °C → add Comins' reagent → vinyl triflate.",
    substrate="Ketone enolates → vinyl triflates. Regiochemistry set by enolization conditions.",
    limitations="Kinetic vs thermodynamic enolate control critical. Triflimide byproduct can complicate purification.",
    mechanism="Enolate O attacks electrophilic S of Tf₂NR → Tf transferred to O → vinyl triflate.",
    refs="Comins, Dehghani, Tetrahedron Lett. 1992, 33, 6299.",
    related="Tf₂O, Heck, Suzuki, Sonogashira",
    functional_groups="ketone, vinyl triflate, enolate")

_nr("Cope Elimination",
    aliases=["Cope elimination"],
    category="Functional Group Interconversion", type_="Elimination",
    summary="Thermal syn-elimination of amine oxides to give alkenes + hydroxylamine.",
    conditions="mCPBA oxidation of amine → amine oxide, then heat (100–150 °C).",
    substrate="Tertiary amines → alkenes. Gives Hofmann (less substituted) product.",
    limitations="High temperature. Requires β-H syn to N-oxide. Hofmann selectivity.",
    mechanism="Concerted syn-periplanar [1,5]-sigmatropic: 5-membered cyclic TS → alkene + R₂NOH.",
    refs="Cope, Foster, JACS 1949, 71, 3929.",
    related="Hofmann Elimination, Chugaev",
    functional_groups="amine, alkene")

_nr("Corey-Chaykovsky Reaction",
    aliases=["Corey-Chaykovsky", "dimethylsulfonium methylide"],
    category="C-C Bond Formation", type_="Cyclopropanation/Epoxidation",
    summary="Dimethylsulfonium methylide or dimethylsulfoxonium methylide → epoxides or cyclopropanes from C=O or C=C.",
    conditions="Me3S⁺I⁻ + NaH (methylide → epoxides). Me3SO⁺I⁻ + NaH (ylide → cyclopropanes from enones).",
    substrate="Aldehydes/ketones + methylide → epoxides. Enones + ylide → cyclopropanes.",
    limitations="Methylide reacts faster (kinetic: epoxides). Ylide stabilized (thermodynamic: cyclopropanes).",
    mechanism="Ylide attacks C=O → betaine → intramolecular displacement → epoxide (or C=C addition → cyclopropane).",
    refs="Corey, Chaykovsky, JACS 1965, 87, 1353.",
    related="Simmons-Smith, Sharpless",
    functional_groups="aldehyde, ketone, epoxide, cyclopropane")

_nr("Curtius Rearrangement",
    aliases=["Curtius rearrangement", "Curtius degradation"],
    category="Rearrangement", type_="Degradation",
    summary="Acyl azides → isocyanates (with loss of N₂) → amines, carbamates, or ureas.",
    conditions="RCON₃ heated (60–100 °C, toluene) or DPPA + Et₃N (one-pot from acid). Trap with ROH → carbamate.",
    substrate="Carboxylic acids → amines (via isocyanate). One carbon shorter (Hofmann-type degradation).",
    limitations="Azide handling hazards. Isocyanate intermediate is moisture-sensitive.",
    mechanism="Acyl azide loses N₂ → nitrene? No: concerted [1,2]-shift (R migrates C→N) → isocyanate.",
    refs="Curtius, Ber. 1890, 23, 3023.",
    related="Hofmann, Lossen, Schmidt rearrangements",
    functional_groups="carboxylic acid, azide, isocyanate, amine")

_nr("Dakin Reaction",
    aliases=["Dakin oxidation"],
    category="Oxidation", type_="Phenol Oxidation",
    summary="H₂O₂ oxidation of o/p-hydroxybenzaldehydes → catechols or hydroquinones.",
    conditions="H₂O₂ (30%), NaOH, H₂O, rt–50 °C.",
    substrate="o-Hydroxybenzaldehyde → catechol. p-Hydroxybenzaldehyde → hydroquinone.",
    limitations="Only works with electron-rich (hydroxylated) arylaldehydes/ketones.",
    mechanism="Baeyer-Villiger-like: HO₂⁻ attacks C=O → Criegee intermediate → aryl migrates → formate ester → hydrolysis → phenol.",
    refs="Dakin, JACS 1909, 42, 477.",
    related="Baeyer-Villiger",
    functional_groups="aldehyde, phenol, catechol")

_nr("Darzens Glycidic Ester Condensation",
    aliases=["Darzens condensation", "Darzens reaction"],
    category="C-C Bond Formation", type_="Condensation",
    summary="α-Halo esters + aldehydes/ketones → glycidic esters (α,β-epoxy esters).",
    conditions="NaOEt or KOtBu, Et₂O or THF, 0 °C → rt.",
    substrate="α-Chloro esters + aldehydes → glycidic esters. Decarboxylation gives aldehydes (homologation).",
    limitations="trans-selective. Competing Claisen condensation with excess base.",
    mechanism="Base → enolate of α-halo ester → aldol with aldehyde → β-hydroxy-α-halo ester → intramolecular SN2 → epoxide.",
    refs="Darzens, Compt. Rend. 1904, 139, 1214.",
    related="Aldol, Reformatsky, Corey-Chaykovsky",
    functional_groups="ester, aldehyde, epoxide")

_nr("Dieckmann Cyclization",
    aliases=["Dieckmann condensation"],
    category="C-C Bond Formation", type_="Intramolecular Condensation",
    summary="Intramolecular Claisen condensation of diesters → cyclic β-keto esters.",
    conditions="NaOEt/EtOH or NaH/THF, reflux. Workup with acid.",
    substrate="Diesters → 5- or 6-membered cyclic β-keto esters. 7-membered rings require high dilution.",
    limitations="5,6-rings only without special conditions. Requires two ester groups.",
    mechanism="Intramolecular enolate attacks second ester → ring closure → elimination of alkoxide.",
    refs="Dieckmann, Ber. 1894, 27, 102.",
    related="Claisen condensation, Aldol",
    functional_groups="ester, β-keto ester, cyclic")

_nr("Doering-LaFlamme Allene Synthesis",
    aliases=["Doering-LaFlamme"],
    category="C-C Bond Formation", type_="Allene Formation",
    summary="Cyclopropanation of alkenes followed by ring opening to allenes using alkyllithium.",
    conditions="CHBr₃/KOtBu → gem-dibromocyclopropane, then nBuLi/Et₂O → allene.",
    substrate="Alkenes → allenes via gem-dihalocyclopropane intermediate.",
    limitations="Requires gem-dihalocyclopropane. Competing elimination side products.",
    mechanism="Carbenoid cyclopropanation → gem-dibromocyclopropane → nBuLi eliminates → carbenoid → ring opens → allene.",
    refs="Doering, LaFlamme, Tetrahedron 1958, 2, 75.",
    related="Simmons-Smith, Corey-Fuchs",
    functional_groups="alkene, allene, cyclopropane")

_nr("Enders SAMP/RAMP Hydrazones",
    aliases=["SAMP hydrazone", "RAMP hydrazone", "Enders asymmetric alkylation"],
    category="C-C Bond Formation", type_="Asymmetric α-Alkylation",
    summary="Chiral hydrazones enable asymmetric α-alkylation of ketones/aldehydes.",
    conditions="Ketone + SAMP → hydrazone, LDA, alkyl halide, then O₃ or MeI/LiAlH₄ to cleave.",
    substrate="Ketones, aldehydes → α-alkylated products with >95% ee.",
    limitations="Requires chiral auxiliary attachment/removal. Multi-step.",
    mechanism="SAMP forms hydrazone → LDA → aza-enolate → alkylation anti to methoxymethyl → cleave auxiliary.",
    refs="Enders, Eichenauer, Angew. Chem. Int. Ed. 1976, 15, 549.",
    related="Evans oxazolidinone, Myers auxiliary",
    functional_groups="ketone, hydrazone, alkyl halide")

_nr("Eschweiler-Clarke Reaction",
    aliases=["Eschweiler-Clarke", "reductive methylation"],
    category="C-N Bond Formation", type_="N-Methylation",
    summary="Exhaustive N-methylation of primary/secondary amines using HCHO/HCO₂H.",
    conditions="HCHO (excess), HCO₂H, 80–100 °C. Gives tertiary amines (no quaternization).",
    substrate="1° amines → 3° amines (dimethylated). 2° amines → 3° amines (monomethylated).",
    limitations="Cannot make quaternary salts. Over-alkylation to NMe₂ with 1° amines.",
    mechanism="Iminium formation with HCHO → formate reduces iminium (Leuckart-Wallach mechanism).",
    refs="Eschweiler, Ber. 1905, 38, 880. Clarke, Gillespie, JACS 1933, 55, 4571.",
    related="Reductive amination, Leuckart-Wallach",
    functional_groups="amine, formaldehyde")

_nr("Evans Aldol Reaction",
    aliases=["Evans aldol", "Evans oxazolidinone"],
    category="C-C Bond Formation", type_="Asymmetric Aldol",
    summary="Chiral oxazolidinone auxiliaries give syn-aldol products with high diastereoselectivity.",
    conditions="N-Acyloxazolidinone + Bu₂BOTf, Et₃N, –78 °C → add RCHO. Then LiOH/H₂O₂ or LiOBn cleavage.",
    substrate="Propionyl oxazolidinone + aldehyde → syn-β-hydroxy acid derivatives (>95:5 dr, >98% ee).",
    limitations="Stoichiometric chiral auxiliary (recoverable). Extra steps for attachment/cleavage.",
    mechanism="B-enolate (Z) → Zimmerman-Traxler 6-membered TS → syn selectivity via chair TS with auxiliary shielding one face.",
    refs="Evans, Bartroli, Shih, JACS 1981, 103, 2127.",
    related="Aldol, Mukaiyama, Crimmins",
    functional_groups="aldehyde, β-hydroxy, oxazolidinone")

_nr("Favorskii Rearrangement",
    aliases=["Favorskii rearrangement"],
    category="Rearrangement", type_="Ring Contraction",
    summary="α-Halo ketones + base → carboxylic acids or esters via cyclopropanone intermediate.",
    conditions="α-Haloketone + NaOH or NaOMe, then acid workup.",
    substrate="α-Haloketones → carboxylic acids (ring contraction if cyclic).",
    limitations="Requires α-halo ketone. Competing elimination.",
    mechanism="Base → enolate → intramolecular displacement of halide → cyclopropanone → nucleophilic ring opening.",
    refs="Favorskii, JRCS 1894, 26, 590.",
    related="Benzilic acid rearrangement",
    functional_groups="ketone, halide, carboxylic acid")

_nr("Ferrier Rearrangement",
    aliases=["Ferrier rearrangement", "Ferrier glycosylation"],
    category="Substitution", type_="Glycosylation",
    summary="Lewis acid-catalyzed allylic rearrangement of glycals → 2,3-unsaturated glycosides.",
    conditions="Glycal + nucleophile + BF₃·Et₂O or SnCl₄, CH₂Cl₂, 0 °C to rt.",
    substrate="Glycals (2,3-unsaturated sugars) + ROH, TMSN₃, etc → C-1 glycosides with C2-C3 unsaturation.",
    limitations="α/β selectivity depends on conditions. Requires glycal substrate.",
    mechanism="Lewis acid activates C-3 OAc leaving group → allylic oxocarbenium → nucleophilic attack at C-1.",
    refs="Ferrier, Prasad, J. Chem. Soc. C 1969, 570.",
    related="Koenigs-Knorr, Schmidt glycosylation",
    functional_groups="glycal, glycoside, alcohol")

_nr("Hantzsch Pyridine Synthesis",
    aliases=["Hantzsch synthesis", "Hantzsch pyridine"],
    category="Heterocycle Formation", type_="Multicomponent",
    summary="One-pot synthesis of 1,4-dihydropyridines from aldehyde + β-ketoester + ammonia.",
    conditions="2 eq β-ketoester + aldehyde + NH₃ (or NH₄OAc), EtOH, reflux. Then oxidize (HNO₃ or DDQ) → pyridine.",
    substrate="Aldehydes + β-keto esters + ammonia → 1,4-dihydropyridines → pyridines.",
    limitations="Limited substitution patterns. Symmetric products from self-condensation.",
    mechanism="Knoevenagel + enamine formation → cyclocondensation → 1,4-DHP → aromatize to pyridine.",
    refs="Hantzsch, Ber. 1881, 14, 1637.",
    related="Chichibabin, Paal-Knorr",
    functional_groups="aldehyde, ester, pyridine")

_nr("Hiyama-Denmark Coupling",
    aliases=["Denmark coupling", "Hiyama-Denmark"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of organosilanes with aryl halides, activated by fluoride or Lewis base.",
    conditions="ArX + R-SiMe₃ or R-Si(OR)₃, Pd(0), TBAF or KF, THF.",
    substrate="Aryl/vinyl halides/triflates + vinyl/aryl silanes → biaryls, dienes.",
    limitations="Requires fluoride activation. Organosilanes less reactive than boronic acids.",
    mechanism="Oxidative addition → fluoride activates Si (pentacoordinate) → transmetalation → reductive elimination.",
    refs="Hiyama, Hatanaka, Pure Appl. Chem. 1994, 66, 1471. Denmark, Regens, Acc. Chem. Res. 2008, 41, 1486.",
    related="Suzuki, Hiyama, Stille",
    functional_groups="aryl halide, silane, biaryl")

_nr("Hofmann Elimination",
    aliases=["Hofmann elimination"],
    category="Functional Group Interconversion", type_="Elimination",
    summary="Quaternary ammonium hydroxides undergo E2 elimination to give less-substituted (Hofmann) alkene.",
    conditions="Amine → MeI (exhaust. methylation) → Me₃N⁺ → Ag₂O/H₂O → heat (100–200 °C).",
    substrate="Quaternary ammonium salts → terminal/less-substituted alkenes (anti-Zaitsev).",
    limitations="Multi-step. High temperature. Competing Zaitsev product.",
    mechanism="E2 with bulky NMe₃ leaving group → anti-periplanar → less-substituted alkene preferred (steric).",
    refs="Hofmann, Annalen 1851, 78, 253.",
    related="Cope Elimination, Zaitsev",
    functional_groups="amine, alkene")

_nr("Horner Reaction",
    aliases=["Horner reaction", "Horner-Wittig"],
    category="C-C Bond Formation", type_="Olefination",
    summary="Phosphine oxide-stabilized carbanions react with aldehydes/ketones → alkenes.",
    conditions="Ph₃P=O-stabilized carbanion + RCHO, nBuLi, THF, –78 °C → rt.",
    substrate="Phosphine oxides + aldehydes → (E)-alkenes. Phosphine oxide byproduct easily removed.",
    limitations="Less reactive than Wittig ylides. Requires strong base.",
    mechanism="Like Wittig but with P(O)Ph₂: betaine → oxaphosphetane → alkene + Ph₂P(O)OH (water soluble).",
    refs="Horner, Hoffmann, Wippel, Ber. 1958, 91, 61.",
    related="HWE, Wittig, Julia Olefination",
    functional_groups="phosphine oxide, aldehyde, alkene")

_nr("Ireland-Claisen Rearrangement",
    aliases=["Ireland-Claisen"],
    category="Rearrangement", type_="[3,3]-Sigmatropic",
    summary="Silyl ketene acetal from ester enolate undergoes [3,3] sigmatropic rearrangement → γ,δ-unsaturated acid.",
    conditions="Ester + LDA → enolate + TMSCl → silyl ketene acetal → heat (60–100 °C) → carboxylic acid.",
    substrate="Allyl esters → γ,δ-unsaturated carboxylic acids. E/Z enolate geometry → syn/anti selectivity.",
    limitations="Requires allyl ester substrate. Stereocontrol via enolate geometry.",
    mechanism="[3,3]-sigmatropic via chair TS: allyl vinyl ether rearrangement. Z-enolate → syn, E-enolate → anti.",
    refs="Ireland, Mueller, Willard, JACS 1976, 98, 2868.",
    related="Claisen, Johnson-Claisen, Eschenmoser-Claisen",
    functional_groups="ester, carboxylic acid, alkene")

_nr("Knoevenagel Condensation",
    aliases=["Knoevenagel"],
    category="C-C Bond Formation", type_="Condensation",
    summary="Aldol-like condensation of aldehydes with active methylene compounds (malonate, cyanoacetate, etc.).",
    conditions="RCHO + CH₂(CO₂Et)₂ or CH₂(CN)CO₂Et, piperidine or pyridine cat., EtOH, reflux. Doebner: decarboxylation in situ.",
    substrate="Aldehydes + active methylene → α,β-unsaturated diesters/nitriles. Doebner variant → cinnamic acids.",
    limitations="Ketones generally unreactive. Needs activated methylene (pKa < 20).",
    mechanism="Base → carbanion → aldol addition → dehydration → alkene. Doebner: in situ decarboxylation → monoester.",
    refs="Knoevenagel, Ber. 1898, 31, 2596.",
    related="Aldol, HWE, Henry",
    functional_groups="aldehyde, malonate, α,β-unsaturated ester")

_nr("Kolbe Electrolysis",
    aliases=["Kolbe electrolysis", "Kolbe reaction"],
    category="C-C Bond Formation", type_="Electrochemistry",
    summary="Electrolytic oxidative decarboxylation of carboxylate salts → alkane dimers (R-R).",
    conditions="RCO₂Na in H₂O or MeOH, Pt electrodes, constant current. Also: Hofer-Moest (forms alcohol/alkene).",
    substrate="Carboxylic acid salts → symmetric alkanes (R-R coupling). Long-chain fatty acids work well.",
    limitations="Symmetric coupling only. Competing Hofer-Moest oxidation.",
    mechanism="Anodic: RCO₂⁻ → RCO₂• → R• + CO₂ → radical coupling → R-R.",
    refs="Kolbe, Annalen 1849, 69, 257.",
    related="Barton decarboxylation",
    functional_groups="carboxylic acid, alkane")

_nr("Kulinkovich Reaction",
    aliases=["Kulinkovich", "Kulinkovich cyclopropanation"],
    category="C-C Bond Formation", type_="Cyclopropanation",
    summary="Ti(OiPr)₄/EtMgBr-mediated conversion of esters → cyclopropanols.",
    conditions="Ester + EtMgBr (2 eq), Ti(OiPr)₄ (cat.), Et₂O, 0 °C → rt.",
    substrate="Esters → 1-substituted cyclopropanols. Also amides (de Meijere variant) → cyclopropylamines.",
    limitations="Only with EtMgBr (or similar). Difficult with sterically hindered esters.",
    mechanism="EtMgBr + Ti(OiPr)₄ → titanacyclopropane → inserts into C=O of ester → oxatitanacycle → ring closure → cyclopropanol.",
    refs="Kulinkovich, Sviridov, Vasilevskii, Synthesis 1991, 234.",
    related="Simmons-Smith, Corey-Chaykovsky",
    functional_groups="ester, cyclopropane, alcohol")

_nr("Lawesson's Reagent",
    aliases=["Lawesson reagent", "thionation"],
    category="Functional Group Interconversion", type_="Thionation",
    summary="C=O → C=S conversion (thionation) of amides, esters, ketones, and lactones.",
    conditions="Lawesson's reagent (0.5–1.0 eq), toluene, 80–110 °C.",
    substrate="Amides → thioamides. Ketones → thioketones. Esters → thioesters. Lactones → thiolactones.",
    limitations="Harsh conditions. Over-thionation possible. Acidic conditions decompose reagent.",
    mechanism="Lawesson's reagent ↔ dithiaphosphetane → [2+2] with C=O → 4-membered ring → cycloreversion → C=S + P=O.",
    refs="Lawesson, Bull. Soc. Chim. Belg. 1978, 87, 229.",
    related="P₄S₁₀",
    functional_groups="amide, ketone, thioamide, thioketone")

_nr("Lossen Rearrangement",
    aliases=["Lossen rearrangement"],
    category="Rearrangement", type_="Degradation",
    summary="Activated hydroxamic acids → isocyanates via [1,2]-shift (like Curtius/Hofmann).",
    conditions="Hydroxamic acid + activating agent (CDI, SOCl₂, p-TsCl) → isocyanate. Trap with ROH/H₂O.",
    substrate="Hydroxamic acids → isocyanates → amines, carbamates, or ureas.",
    limitations="Requires pre-formed hydroxamic acid + activation.",
    mechanism="Activation of O-H → departure of leaving group → [1,2]-shift of R from C to N → isocyanate + byproduct.",
    refs="Lossen, Annalen 1872, 161, 347.",
    related="Curtius, Hofmann, Schmidt",
    functional_groups="hydroxamic acid, isocyanate, amine")

_nr("Luche Reduction",
    aliases=["Luche reduction"],
    category="Reduction", type_="Selective Reduction",
    summary="1,2-selective reduction of α,β-unsaturated ketones using NaBH₄/CeCl₃.",
    conditions="NaBH₄, CeCl₃·7H₂O, MeOH, 0 °C → rt. Quick reaction (minutes).",
    substrate="Enones → allylic alcohols (1,2-reduction). Selective over 1,4-reduction.",
    limitations="Only for α,β-unsaturated carbonyls. Not for isolated C=C.",
    mechanism="CeCl₃ activates C=O (hard Lewis acid) → NaBH₄ delivers H⁻ to C=O (1,2) instead of conjugate (1,4).",
    refs="Luche, JACS 1978, 100, 2226.",
    related="NaBH₄, CBS, Meerwein-Ponndorf-Verley",
    functional_groups="enone, allylic alcohol")

_nr("Matteson Homologation",
    aliases=["Matteson homologation"],
    category="Homologation", type_="Boron Homologation",
    summary="Iterative one-carbon homologation of boronic esters with high stereocontrol.",
    conditions="R-B(pin) + CHCl₂Li (or ICH₂Cl/nBuLi) → α-chloroboronate → ZnCl₂/MgBr₂ stereospecific substitution.",
    substrate="Boronic esters → homologated boronic esters, iteratively with stereocontrol.",
    limitations="Requires careful handling of α-halocarbanions. Multi-step for complex targets.",
    mechanism="Carbenoid inserts into C-B bond → α-haloboronate → 1,2-migration with inversion → new C-C bond.",
    refs="Matteson, Mah, JACS 1982, 104, 4471.",
    related="Suzuki, Aggarwal homologation",
    functional_groups="boronic ester, alkyl halide")

_nr("Michael Reaction",
    aliases=["Michael addition", "conjugate addition", "1,4-addition"],
    category="C-C Bond Formation", type_="Conjugate Addition",
    summary="Nucleophilic 1,4-addition to α,β-unsaturated carbonyl compounds.",
    conditions="Nucleophile + enone/enal, base (NaOEt, DBU, or organocatalyst), or CuLn for cuprate 1,4-addition.",
    substrate="Malonates, enolates, or cuprates + enones → 1,5-dicarbonyls. Also aza-Michael (amines), oxa-Michael (alcohols).",
    limitations="Competing 1,2-addition with hard nucleophiles. Regioselectivity control needed.",
    mechanism="Deprotonation → carbanion attacks β-carbon of enone → enolate intermediate → protonation at α.",
    refs="Michael, J. Prakt. Chem. 1887, 35, 349.",
    related="Robinson annulation, Stetter, Mukaiyama-Michael",
    functional_groups="enone, malonate, 1,5-dicarbonyl")

_nr("Minisci Reaction",
    aliases=["Minisci reaction", "Minisci radical addition"],
    category="Radical", type_="Radical C-H Functionalization",
    summary="Radical addition to protonated heteroaromatics (pyridine, quinoline, etc.) for C-H functionalization.",
    conditions="Hetarene + RCO₃H or R-I + AgNO₃/persulfate, H₂SO₄ (acid), 60–80 °C.",
    substrate="Pyridines, quinolines, isoquinolines → C-2/C-4 alkylated/acylated products. C-2 preferred.",
    limitations="Regioselectivity can be poor. Mixtures of mono/disubstitution. Requires acid medium.",
    mechanism="Radical generated (from acid, peroxide, or photoredox) → adds to protonated heterocycle → rearomatization.",
    refs="Minisci, Bernardi, Bertini, Synthesis 1971, 650.",
    related="Giese radical addition",
    functional_groups="pyridine, radical, carboxylic acid")

_nr("Mukaiyama Aldol Reaction",
    aliases=["Mukaiyama aldol", "Mukaiyama"],
    category="C-C Bond Formation", type_="Lewis Acid Aldol",
    summary="Lewis acid-catalyzed aldol of silyl enol ethers with aldehydes.",
    conditions="Silyl enol ether + RCHO + TiCl₄ or BF₃·Et₂O, CH₂Cl₂, –78 °C. Asymmetric: chiral Cu(II), Sn(II), or B catalysts.",
    substrate="Silyl enol ethers + aldehydes → β-hydroxy ketones. Vinylogous variant for γ-addition.",
    limitations="Requires pre-formed silyl enol ether. Lewis acid can racemize product.",
    mechanism="Lewis acid activates aldehyde → nucleophilic addition of silyl enol ether → β-silyloxy ketone → deprotect.",
    refs="Mukaiyama, Narasaka, Banno, Chem. Lett. 1973, 2, 1011.",
    related="Aldol, Evans, Masamune",
    functional_groups="silyl enol ether, aldehyde, β-hydroxy ketone")

_nr("Myers Asymmetric Alkylation",
    aliases=["Myers alkylation", "pseudoephedrine auxiliary"],
    category="C-C Bond Formation", type_="Asymmetric Alkylation",
    summary="Pseudoephedrine amide auxiliaries for highly enantioselective α-alkylation.",
    conditions="Pseudoephedrine amide + LDA + LiCl, THF, 0 °C → alkyl halide. Cleave: LiOH/H₂O₂ or LiNH₂/NH₃(l).",
    substrate="Amides → α-alkylated carboxylic acids with >95% ee. Broad substrate scope.",
    limitations="Stoichiometric auxiliary. Extra steps for attachment/cleavage.",
    mechanism="LDA/LiCl → chelated (Z)-enolate → alkylation on exposed face → high ee.",
    refs="Myers, Yang, Chen, Cebulak, JACS 1997, 119, 6496.",
    related="Evans, Enders, Oppolzer",
    functional_groups="amide, alkyl halide, carboxylic acid")

_nr("Negishi Coupling",
    aliases=["Negishi coupling"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd or Ni-catalyzed cross-coupling of organozinc reagents with organic halides.",
    conditions="R-ZnX + R'-X, Pd(PPh₃)₄ or Pd(dba)₂/ligand, THF, 0 °C → rt. Very high functional group tolerance.",
    substrate="Aryl, vinyl, alkyl zinc + aryl/vinyl halides/triflates → coupled products.",
    limitations="Organozinc preparation requires anhydrous conditions. Less atom-economical than Suzuki.",
    mechanism="Oxidative addition of R'-X to Pd(0) → transmetalation with R-ZnX → reductive elimination → R-R' + Pd(0).",
    refs="Negishi, King, Okukado, J. Org. Chem. 1977, 42, 1821.",
    related="Suzuki, Kumada, Stille, Heck",
    functional_groups="organozinc, aryl halide, biaryl")

_nr("Nozaki-Hiyama-Kishi Reaction",
    aliases=["NHK reaction", "Nozaki-Hiyama-Kishi"],
    category="C-C Bond Formation", type_="Allylation/Carbonyl Addition",
    summary="CrCl₂-mediated coupling of vinyl/allyl halides with aldehydes (highly chemoselective for aldehyde over ketone).",
    conditions="CrCl₂ (excess, 2–6 eq) + NiCl₂ (cat.), DMF, rt. Modern: Fürstner's cat. Cr + Mn reductant.",
    substrate="Vinyl iodides/bromides + aldehydes → allylic alcohols. Remarkable aldehyde selectivity.",
    limitations="Stoichiometric Cr (toxic, environmentally problematic). Catalytic versions exist but less developed.",
    mechanism="Cr(II) reduces vinyl halide → vinyl-Cr(III) → nucleophilic addition to aldehyde → allylic alcohol + Cr(III) salt.",
    refs="Nozaki, Oshima, Takai, JACS 1977, 99, 3532. Kishi: application to palytoxin.",
    related="Grignard, Barbier",
    functional_groups="vinyl halide, aldehyde, allylic alcohol")

_nr("Passerini Reaction",
    aliases=["Passerini reaction", "Passerini 3CR"],
    category="C-C Bond Formation", type_="Multicomponent",
    summary="Three-component reaction of isocyanide + aldehyde + carboxylic acid → α-acyloxyamide.",
    conditions="R-NC + R'CHO + R''CO₂H, CH₂Cl₂ or MeOH, rt, 24–48 h.",
    substrate="Isocyanides + aldehydes + carboxylic acids → α-acyloxyamides. Huge diversity-oriented potential.",
    limitations="Isocyanide odor. Limited stereocontrol (racemic).",
    mechanism="Aldehyde + acid → H-bonded complex → isocyanide inserts into C=O → α-adduct → acyl migration → product.",
    refs="Passerini, Gazz. Chim. Ital. 1921, 51, 126.",
    related="Ugi, Biginelli",
    functional_groups="isocyanide, aldehyde, carboxylic acid, amide")

_nr("Paternò-Büchi Reaction",
    aliases=["Paterno-Buchi", "Paternò-Büchi"],
    category="Cycloaddition", type_="[2+2] Photocycloaddition",
    summary="Photochemical [2+2] cycloaddition of carbonyl compounds with alkenes → oxetanes.",
    conditions="C=O + alkene, hν (UV, 254–350 nm), solvent (benzene, acetone, CH₂Cl₂).",
    substrate="Aldehydes/ketones + electron-rich alkenes → oxetanes. Regiochemistry: head-to-head preferred.",
    limitations="Requires UV irradiation. Competing [2+2] at C=C. Product strain (4-membered ring).",
    mechanism="Carbonyl excited (S₁ or T₁ via ISC) → addition to alkene → diradical → ring closure → oxetane.",
    refs="Paternò, Chieffi, Gazz. Chim. Ital. 1909, 39, 341. Büchi, Inman, Lipinsky, JACS 1954, 76, 4327.",
    related="[2+2] photocycloaddition, Norrish",
    functional_groups="carbonyl, alkene, oxetane")

_nr("Petasis Reaction",
    aliases=["Petasis borono-Mannich", "Petasis reaction"],
    category="C-C Bond Formation", type_="Multicomponent",
    summary="Three-component coupling of amine + aldehyde + boronic acid → substituted amine.",
    conditions="R₂NH + RCHO + ArB(OH)₂, CH₂Cl₂ or MeOH, rt. No metal catalyst needed.",
    substrate="Amines + glyoxylic acid or salicylaldehydes + boronic acids → α-amino acids, aminophenols.",
    limitations="Best with glyoxylic acid or salicylaldehydes (need α-hydroxy). Aliphatic boronic acids less reactive.",
    mechanism="Imine/iminium formation → boronate coordinates to adjacent OH → intramolecular transfer of aryl from B to C.",
    refs="Petasis, Akritopoulou, Tetrahedron Lett. 1993, 34, 583.",
    related="Mannich, Ugi, Suzuki",
    functional_groups="amine, aldehyde, boronic acid")

_nr("Pictet-Spengler Reaction",
    aliases=["Pictet-Spengler"],
    category="Heterocycle Formation", type_="Ring Closure",
    summary="Acid-catalyzed condensation of β-arylethylamines with aldehydes → tetrahydroisoquinolines/β-carbolines.",
    conditions="Tryptamine or PEA derivative + RCHO, TFA or AcOH, CH₂Cl₂ or MeOH, rt–reflux.",
    substrate="Tryptamines + aldehydes → β-carbolines. PEA + aldehydes → tetrahydroisoquinolines.",
    limitations="Electron-poor arenes are sluggish. Stereocontrol challenging (chiral phosphoric acid catalysts).",
    mechanism="Iminium formation → electrophilic cyclization → Wheland intermediate → rearomatization.",
    refs="Pictet, Spengler, Ber. 1911, 44, 2030.",
    related="Bischler-Napieralski, Mannich",
    functional_groups="amine, aldehyde, isoquinoline, β-carboline")

_nr("Pinner Reaction",
    aliases=["Pinner synthesis", "Pinner reaction"],
    category="Functional Group Interconversion", type_="Nitrile Conversion",
    summary="Acid-catalyzed conversion of nitriles to imino esters (Pinner salts), then to esters, amidines, or amides.",
    conditions="R-CN + ROH + dry HCl(g) → Pinner salt (imidate·HCl). Then: H₂O → ester, NH₃ → amidine, R'OH → orthoester.",
    substrate="Nitriles → imino esters → esters, amidines, amides, orthoesters.",
    limitations="Anhydrous HCl required. Pinner salt hydrolysis must be controlled.",
    mechanism="HCl protonates nitrile → alcohol attacks activated N≡C → iminium → Pinner salt. Hydrolysis → ester.",
    refs="Pinner, Klein, Ber. 1877, 10, 1889.",
    related="Ritter reaction",
    functional_groups="nitrile, ester, amidine")

_nr("Prins Reaction",
    aliases=["Prins reaction", "Prins cyclization"],
    category="C-C Bond Formation", type_="Carbonyl-Ene/Cyclization",
    summary="Acid-catalyzed addition of aldehydes to alkenes → homoallylic alcohols, 1,3-dioxanes, or tetrahydropyrans.",
    conditions="R-CH=CH₂ + RCHO, Lewis acid (BF₃, SnCl₄, InBr₃), CH₂Cl₂ or neat, –78 °C to rt.",
    substrate="Alkenes + aldehydes → homoallylic alcohols (kinetic) or tetrahydropyrans (cyclization). Prins-type cascade for natural products.",
    limitations="Product depends on conditions (alcohol vs THP vs dioxane). Competing polymerization.",
    mechanism="Lewis acid activates aldehyde → ene reaction or carbocation formation → C-C bond → THP formation or trapping.",
    refs="Prins, Chem. Weekblad 1919, 16, 1072.",
    related="Conia-ene, Mukaiyama",
    functional_groups="alkene, aldehyde, tetrahydropyran, alcohol")

_nr("Ramberg-Bäcklund Reaction",
    aliases=["Ramberg-Backlund", "Ramberg-Bäcklund"],
    category="C-C Bond Formation", type_="Sulfone Elimination",
    summary="α-Halo sulfones → alkenes via episulfone intermediate with extrusion of SO₂.",
    conditions="α-Haloalkylsulfone + KOtBu or NaOH, or Meyers' conditions (KOH/Al₂O₃, CHCl₃, rt).",
    substrate="α-Haloalkylsulfones → alkenes (with loss of SO₂). Meyers conditions very mild.",
    limitations="Requires α-halo sulfone precursor. E/Z control limited.",
    mechanism="Base abstracts α-H → carbanion displaces halide → episulfone → cheletropic loss of SO₂ → alkene.",
    refs="Ramberg, Bäcklund, Arkiv Kemi 1940, 13A, 1.",
    related="Julia Olefination, Peterson",
    functional_groups="sulfone, alkene, halide")

_nr("Ring-Closing Metathesis",
    aliases=["RCM", "ring-closing metathesis"],
    category="C-C Bond Formation", type_="Metathesis",
    summary="Ru-catalyzed intramolecular olefin metathesis to form carbocyclic/heterocyclic rings.",
    conditions="Grubbs I (PCy₃), Grubbs II (NHC), or Hoveyda-Grubbs cat., CH₂Cl₂ or toluene, 40 °C, dilute. Ethylene removed under vacuum/N₂.",
    substrate="Dienes → 5–8-membered rings (best for 5, 6). Macrocycles also possible at high dilution.",
    limitations="Ring strain limits small rings (<5). Catalyst poisoning by N, S, P. Dilution needed for macro-RCM.",
    mechanism="Ru=CHR → [2+2] with terminal alkene → metallacyclobutane → retro [2+2] → new Ru=CH → intramolecular closure → product + CH₂=CH₂.",
    refs="Grubbs, JACS 1995, 117, 5503. Schrock, Murdzek, JACS 1990, 112, 3875.",
    related="Olefin metathesis, ROMP, CM",
    functional_groups="diene, alkene, ring")

_nr("Ritter Reaction",
    aliases=["Ritter reaction"],
    category="C-N Bond Formation", type_="Amide Formation",
    summary="Acid-catalyzed addition of nitriles to carbocations (from alkenes/alcohols) → amides.",
    conditions="R-OH or R-CH=CH₂ + R'CN, H₂SO₄ or BF₃·Et₂O, 0 °C → rt. Then H₂O workup → amide.",
    substrate="Tertiary alcohols/alkenes + nitriles → N-tert-alkyl amides. Good for tert-amines (hydrolysis).",
    limitations="Requires stable carbocation (3°, benzylic, allylic). 1° and 2° give poor results.",
    mechanism="Acid → carbocation → attacks nitrile N → nitrilium → H₂O adds → amide (tautomerize).",
    refs="Ritter, Minieri, JACS 1948, 70, 4045.",
    related="Pinner, Beckmann",
    functional_groups="nitrile, alcohol, amide")

_nr("Robinson Annulation",
    aliases=["Robinson annulation"],
    category="C-C Bond Formation", type_="Annulation",
    summary="Tandem Michael addition + intramolecular aldol → cyclohexenone ring construction.",
    conditions="Ketone + methyl vinyl ketone, NaOH or KOH, EtOH/H₂O, rt → reflux. Also: Hajos-Parrish (proline-catalyzed).",
    substrate="Ketones + MVK → 2-cyclohexenones. Key step in steroid synthesis (Hajos-Parrish-Eder-Sauer-Wiechert).",
    limitations="Polymerization of MVK. Regiochemistry with unsymmetrical ketones.",
    mechanism="Michael addition of enolate to enone → 1,5-diketone → intramolecular aldol → β-hydroxy ketone → dehydration → cyclohexenone.",
    refs="Robinson, J. Chem. Soc. 1935, 1604.",
    related="Michael, Aldol, Hajos-Parrish",
    functional_groups="ketone, enone, cyclohexenone")

_nr("Roush Allylboration",
    aliases=["Roush allylboration", "allylboration"],
    category="C-C Bond Formation", type_="Allylation",
    summary="Chiral allylboronates add to aldehydes with high diastereo- and enantioselectivity.",
    conditions="Chiral tartrate allylboronate + RCHO, toluene, –78 °C. Also: Brown (Ipc₂B), Roush (DIPT tartrate).",
    substrate="Aldehydes → homoallylic alcohols with up to >95% ee. E-crotyl → anti, Z-crotyl → syn.",
    limitations="Requires chiral auxiliary on boron. Stoichiometric chiral reagent.",
    mechanism="6-membered Zimmerman-Traxler TS: aldehyde coordinates B, allyl transferred via chair-like TS → homoallylic alcohol.",
    refs="Roush, Walts, Hoong, JACS 1985, 107, 8186.",
    related="Brown allylation, Matteson, crotylation",
    functional_groups="aldehyde, allylboron, homoallylic alcohol")

_nr("Sakurai-Hosomi Allylation",
    aliases=["Sakurai allylation", "Hosomi-Sakurai"],
    category="C-C Bond Formation", type_="Allylation",
    summary="Lewis acid-mediated allylation of aldehydes/ketones with allylsilanes.",
    conditions="Allyl-TMS + RCHO + TiCl₄ or BF₃·Et₂O, CH₂Cl₂, –78 °C.",
    substrate="Aldehydes, ketones, acetals → homoallylic alcohols. Good functional group tolerance.",
    limitations="Requires Lewis acid. Less enantioselective than allylborations (unless chiral Lewis acid).",
    mechanism="Lewis acid activates C=O → allylsilane attacks → SE2' → homoallylic alcohol + TMS-X.",
    refs="Hosomi, Sakurai, Tetrahedron Lett. 1976, 17, 1295.",
    related="Crotylation, Roush, Mukaiyama aldol",
    functional_groups="allylsilane, aldehyde, homoallylic alcohol")

_nr("Simmons-Smith Cyclopropanation",
    aliases=["Simmons-Smith", "Simmons-Smith reaction"],
    category="Cycloaddition", type_="Cyclopropanation",
    summary="Zn-carbenoid cyclopropanation of alkenes using CH₂I₂/Zn-Cu.",
    conditions="CH₂I₂ + Zn-Cu couple or Et₂Zn (Furukawa mod.), Et₂O or CH₂Cl₂, 0 °C → rt. Charette: directed with chiral bisoxazoline.",
    substrate="Alkenes → cyclopropanes. Directed by adjacent OH (syn delivery). Stereospecific (retention of alkene geometry).",
    limitations="Methylene transfer only (no substituted carbenoids easily). Zn waste.",
    mechanism="ICH₂ZnI (Zn carbenoid) → [2+1] delivery of CH₂ to alkene face → cyclopropane. Butterfly/concerted TS.",
    refs="Simmons, Smith, JACS 1958, 80, 5323.",
    related="Corey-Chaykovsky, Kulinkovich, Doering-LaFlamme",
    functional_groups="alkene, cyclopropane")

_nr("Stetter Reaction",
    aliases=["Stetter reaction", "NHC-catalyzed 1,4-addition"],
    category="C-C Bond Formation", type_="Umpolung Conjugate Addition",
    summary="NHC-catalyzed umpolung: aldehyde → acyl anion equivalent → 1,4-addition to Michael acceptors → 1,4-diketones.",
    conditions="Aldehyde + enone + thiazolium/triazolium salt (NHC) + Et₃N, EtOH, reflux. Intramolecular: Stetter cyclization.",
    substrate="Aldehydes + Michael acceptors → 1,4-dicarbonyls. Intramolecular for chromanones, etc.",
    limitations="Intermolecular often low yielding. Self-condensation of aldehyde (benzoin). Asymmetric variants developing.",
    mechanism="NHC adds to aldehyde → Breslow intermediate (acyl anion equiv.) → conjugate addition → 1,4-product → NHC released.",
    refs="Stetter, Angew. Chem. Int. Ed. 1976, 15, 639.",
    related="Benzoin condensation, Michael",
    functional_groups="aldehyde, enone, 1,4-diketone")

_nr("Stork Enamine Alkylation",
    aliases=["Stork enamine", "enamine alkylation"],
    category="C-C Bond Formation", type_="α-Alkylation",
    summary="Enamines as nucleophiles for α-alkylation/acylation of ketones under mild conditions.",
    conditions="Ketone + 2° amine (pyrrolidine) → enamine → alkyl/acyl halide → hydrolysis → α-substituted ketone.",
    substrate="Ketones → α-alkylated or α-acylated ketones. Regioselective (less-substituted enamine preferred).",
    limitations="Less-substituted enamine forms preferentially (kinetic). Competing N-alkylation.",
    mechanism="Amine + ketone → enamine (C=C-N) → C-alkylation with electrophile → iminium → hydrolysis → ketone.",
    refs="Stork, Terrell, Szmuszkovicz, JACS 1954, 76, 2029.",
    related="Enders SAMP/RAMP, organocatalysis",
    functional_groups="ketone, amine, alkyl halide")

_nr("Strecker Synthesis",
    aliases=["Strecker synthesis", "Strecker amino acid synthesis"],
    category="C-C Bond Formation", type_="Multicomponent",
    summary="Aldehyde + ammonia + HCN → α-aminonitrile → α-amino acid (hydrolysis).",
    conditions="RCHO + NH₄Cl + NaCN (or TMSCN), H₂O or MeOH, 0 °C → rt. Asymmetric: chiral thiourea or phosphoric acid cat.",
    substrate="Aldehydes → α-amino acids (after hydrolysis of nitrile). Asymmetric variants give high ee.",
    limitations="HCN toxicity. Racemic unless chiral catalyst used.",
    mechanism="Imine formation (RCHO + NH₃) → cyanide attacks imine → α-aminonitrile → acid hydrolysis → amino acid.",
    refs="Strecker, Annalen 1850, 75, 27.",
    related="Ugi, Passerini",
    functional_groups="aldehyde, amine, amino acid, nitrile")

_nr("Takai Olefination",
    aliases=["Takai olefination", "Takai reaction"],
    category="C-C Bond Formation", type_="Olefination",
    summary="CrCl₂/CHI₃ converts aldehydes to (E)-vinyl iodides selectively.",
    conditions="RCHO + CHI₃, CrCl₂ (excess), THF, 0 °C → rt. Also CHBr₃ for vinyl bromides.",
    substrate="Aldehydes → (E)-vinyl iodides. Excellent E-selectivity (>95:5).",
    limitations="Stoichiometric Cr(II) (toxic). Not for ketones.",
    mechanism="CrCl₂ reduces CHI₃ → Cr-carbenoid → adds to aldehyde → vinyl metal → β-elimination → (E)-vinyl iodide.",
    refs="Takai, Nitta, Utimoto, JACS 1986, 108, 7408.",
    related="NHK, Wittig, Still-Gennari",
    functional_groups="aldehyde, vinyl iodide")

_nr("Tiffeneau-Demjanov Rearrangement",
    aliases=["Tiffeneau-Demjanov", "ring expansion"],
    category="Rearrangement", type_="Ring Expansion",
    summary="Cyclic α-amino ketones → ring-expanded ketones (+1 carbon) via diazotization.",
    conditions="α-Amino ketone + NaNO₂/AcOH → α-diazo ketone → N₂ loss → ring expansion.",
    substrate="Cyclopentanone → cyclohexanone, cyclohexanone → cycloheptanone, etc.",
    limitations="Competing Wolff rearrangement. Requires α-amino ketone.",
    mechanism="Diazotization → α-diazo ketone → N₂ leaves → 1,2-shift of C-C bond adjacent to carbocation → ring-expanded ketone.",
    refs="Tiffeneau, Weill, Tchoubar, Compt. Rend. 1937, 205, 54.",
    related="Beckmann, Baeyer-Villiger, Buchner",
    functional_groups="ketone, amine, ring expansion")

_nr("Ullmann Coupling",
    aliases=["Ullmann coupling", "Ullmann reaction"],
    category="C-C Bond Formation", type_="Cu-Mediated Coupling",
    summary="Cu-mediated homo- or cross-coupling of aryl halides. Modern: Cu-catalyzed C-N, C-O, C-S coupling (Ullmann-type).",
    conditions="Classical: 2 ArX + Cu powder, 200 °C. Modern: ArX + NuH, CuI (10%), diamine ligand, Cs₂CO₃, DMSO, 80 °C.",
    substrate="Aryl halides → biaryls (classical). ArX + amines/phenols/thiols → C-N/C-O/C-S (modern Ullmann-Goldberg-Buchwald).",
    limitations="Classical requires high T. Modern requires good ligand design for challenging substrates.",
    mechanism="Classical: Cu(0) → oxidative addition → Cu(III)-Ar → second ArX → reductive elimination → Ar-Ar. Modern: Cu(I)/Cu(III) cycle.",
    refs="Ullmann, Bielecki, Ber. 1901, 34, 2174.",
    related="Buchwald-Hartwig, Chan-Lam, Suzuki",
    functional_groups="aryl halide, biaryl, amine, phenol")

_nr("Ugi Reaction",
    aliases=["Ugi reaction", "Ugi 4CR", "Ugi multicomponent"],
    category="C-C Bond Formation", type_="Multicomponent",
    summary="Four-component reaction: amine + aldehyde + carboxylic acid + isocyanide → α-acylaminoamide.",
    conditions="R₂NH + R'CHO + R''CO₂H + R'''NC, MeOH, rt, 24–48 h. Also: Ugi-3CR, Ugi-azide, Ugi-Heck cascades.",
    substrate="Generates enormous molecular diversity from simple building blocks. Peptoid backbones, macrocycles.",
    limitations="Isocyanide odor/toxicity. Racemic (chiral variants under development). Purification challenging.",
    mechanism="Imine formation → protonation by acid → isocyanide attacks iminium → α-adduct → Mumm rearrangement → product.",
    refs="Ugi, Steinbrückner, Angew. Chem. 1960, 72, 267.",
    related="Passerini, Biginelli, Strecker",
    functional_groups="amine, aldehyde, carboxylic acid, isocyanide, amide")

_nr("Vilsmeier-Haack Formylation",
    aliases=["Vilsmeier-Haack", "Vilsmeier formylation"],
    category="C-C Bond Formation", type_="Electrophilic Formylation",
    summary="POCl₃/DMF-mediated formylation of electron-rich arenes (phenols, anilines, heterocycles).",
    conditions="DMF + POCl₃ (0 °C, form Vilsmeier reagent) → add substrate, 60–80 °C → aqueous workup → aldehyde.",
    substrate="Electron-rich aromatics → aryl aldehydes. Pyrroles, indoles, phenol ethers, anilines.",
    limitations="Fails with electron-poor arenes. Free phenols can give Duff or Reimer-Tiemann products.",
    mechanism="POCl₃ + DMF → chloroiminium salt (Vilsmeier reagent) → electrophilic substitution at electron-rich C → iminium → hydrolysis → aldehyde.",
    refs="Vilsmeier, Haack, Ber. 1927, 60, 119.",
    related="Gattermann, Duff, Reimer-Tiemann",
    functional_groups="arene, aldehyde, DMF")

_nr("Weinreb Amide",
    aliases=["Weinreb amide", "Weinreb ketone synthesis"],
    category="C-C Bond Formation", type_="Ketone Synthesis",
    summary="N-Methoxy-N-methylamides (Weinreb amides) react with Grignard/organolithium → ketones (no over-addition).",
    conditions="R-C(=O)N(OMe)Me + R'MgBr (or R'Li), THF, –78 °C → 0 °C. Single equivalent of nucleophile.",
    substrate="Weinreb amides + RMgBr → ketones. Also: + DIBAL → aldehydes.",
    limitations="Requires Weinreb amide preparation (extra step). Sensitive to steric hindrance.",
    mechanism="Nucleophile attacks C=O → tetrahedral intermediate chelated by N-OMe → stable → does not collapse until workup → ketone.",
    refs="Nahm, Weinreb, Tetrahedron Lett. 1981, 22, 3815.",
    related="Grignard + acid chloride, Fukuyama coupling",
    functional_groups="amide, Grignard, ketone, aldehyde")

_nr("Wharton Fragmentation",
    aliases=["Wharton reaction", "Wharton transposition"],
    category="Rearrangement", type_="Transposition",
    summary="α,β-Epoxyketones + hydrazine → allylic alcohols via [2,3]-sigmatropic shift.",
    conditions="α,β-Epoxyketone + H₂NNH₂·H₂O, EtOH/AcOH, reflux.",
    substrate="α,β-Epoxyketones → allylic alcohols (transposition of C=O).",
    limitations="Requires α,β-epoxy ketone substrate specifically.",
    mechanism="Hydrazone forms → [2,3]-sigmatropic shift → diazene intermediate → loss of N₂ → allylic alcohol.",
    refs="Wharton, Bohlen, J. Org. Chem. 1961, 26, 3615.",
    related="Bamford-Stevens, Shapiro",
    functional_groups="epoxide, ketone, allylic alcohol")

_nr("Wittig Reaction",
    aliases=["Wittig reaction", "Wittig olefination"],
    category="C-C Bond Formation", type_="Olefination",
    summary="Phosphonium ylide + aldehyde/ketone → alkene + Ph₃P=O.",
    conditions="Ph₃P=CHR (ylide from Ph₃P⁺CH₂R X⁻ + base) + RCHO, THF, –78 °C → rt. Non-stabilized → Z, stabilized → E.",
    substrate="Aldehydes + ylides → alkenes. Non-stabilized (R = alkyl) → cis. Stabilized (R = CO₂R, CN) → trans.",
    limitations="Triphenylphosphine oxide removal. Z/E selectivity depends on ylide type.",
    mechanism="Ylide + aldehyde → betaine → oxaphosphetane (via [2+2] or stepwise) → retro-[2+2] → alkene + Ph₃P=O.",
    refs="Wittig, Geissler, Annalen 1953, 580, 44.",
    related="HWE, Still-Gennari, Julia, Tebbe",
    functional_groups="aldehyde, ylide, alkene")

_nr("Wolff Rearrangement",
    aliases=["Wolff rearrangement", "Arndt-Eistert (Wolff step)"],
    category="Rearrangement", type_="Carbene Rearrangement",
    summary="α-Diazo ketones → ketenes via Wolff rearrangement (loss of N₂, [1,2]-shift).",
    conditions="α-Diazoketone + hν or Ag₂O or heat → ketene → trap with H₂O, ROH, or amine.",
    substrate="Acyl chloride + CH₂N₂ → α-diazoketone → ketene → acid, ester, or amide (one carbon homologation).",
    limitations="Diazoketone hazards. Wolff competes with C-H insertion (Rh catalysis).",
    mechanism="α-Diazoketone → loss of N₂ → carbene/carbenoid → [1,2]-shift → ketene → nucleophilic trapping.",
    refs="Wolff, Annalen 1912, 394, 23.",
    related="Arndt-Eistert, Curtius",
    functional_groups="diazo, ketene, carboxylic acid")

_nr("Yamaguchi Esterification",
    aliases=["Yamaguchi esterification", "Yamaguchi macrolactonization"],
    category="Functional Group Interconversion", type_="Esterification/Macrolactonization",
    summary="2,4,6-Trichlorobenzoyl chloride-mediated esterification, especially for macrolactonization.",
    conditions="RCO₂H + 2,4,6-Cl₃C₆H₂COCl, Et₃N, THF → mixed anhydride → add ROH + DMAP → ester. Macrolactonization: high dilution.",
    substrate="Seco-acids → macrolactones (12–20-membered). Also intermolecular esterification of hindered alcohols.",
    limitations="High dilution for macrolactonization. Expensive reagent. DMAP crucial for rate.",
    mechanism="Acid + trichlorobenzoyl chloride → mixed anhydride → DMAP-catalyzed aminolysis → acyl-DMAP → alcoholysis → ester.",
    refs="Inanaga, Hirata, Saeki, Katsuki, Yamaguchi, Bull. Chem. Soc. Jpn. 1979, 52, 1989.",
    related="Steglich, Fischer, Mitsunobu, Shiina",
    functional_groups="carboxylic acid, alcohol, ester, macrolactone")

_nr("Zincke Aldehyde",
    aliases=["Zincke aldehyde", "Zincke reaction"],
    category="Heterocycle Formation", type_="Ring Opening/Closing",
    summary="Pyridinium salts (Zincke salts) react with primary amines → N-substituted pyridinium salts (pyridine N-exchange).",
    conditions="Pyridine + 2,4-dinitro-chlorobenzene → Zincke salt → R-NH₂ → N-substituted pyridinium → base → Zincke aldehyde (5-amino-2,4-pentadienal).",
    substrate="Pyridinium salts + 1° amines → pyridinium exchange. Ring-opened intermediates = Zincke aldehydes.",
    limitations="Requires DNB activation. Multi-step. Low yields for electron-poor amines.",
    mechanism="Ring opens via retro-electrocyclic → pentadienal-iminium (Zincke aldehyde) → new amine condenses → ring closes → new pyridinium.",
    refs="Zincke, Annalen 1903, 330, 361.",
    related="Chichibabin",
    functional_groups="pyridine, amine")

_nr("Biginelli Reaction",
    aliases=["Biginelli reaction", "Biginelli 3CR"],
    category="Heterocycle Formation", type_="Multicomponent",
    summary="Three-component synthesis of dihydropyrimidinones (DHPMs) from aldehyde + urea + β-ketoester.",
    conditions="RCHO + urea + β-ketoester, HCl or Lewis acid (Yb(OTf)₃, FeCl₃), EtOH, reflux.",
    substrate="Aldehydes + urea + 1,3-dicarbonyls → 3,4-dihydropyrimidin-2(1H)-ones. Enantioselective with chiral phosphoric acid.",
    limitations="Moderate yields with aliphatic aldehydes. Slow with bulky substrates.",
    mechanism="Aldehyde + urea → N-acyliminium → β-ketoester → Mannich-type → cyclization → DHPM.",
    refs="Biginelli, Ber. 1891, 24, 1317.",
    related="Hantzsch, Ugi, Passerini",
    functional_groups="aldehyde, urea, β-ketoester, dihydropyrimidinone")

_nr("Corey-Bakshi-Shibata Reduction",
    aliases=["CBS reduction", "Corey-Bakshi-Shibata"],
    category="Reduction", type_="Asymmetric Reduction",
    summary="Oxazaborolidine-catalyzed asymmetric reduction of prochiral ketones with BH₃.",
    conditions="Ketone + BH₃·THF (or catecholborane), CBS catalyst (10 mol%), toluene, –20 °C to rt.",
    substrate="Prochiral ketones → chiral secondary alcohols with >95% ee. Aryl-alkyl ketones best.",
    limitations="Both enantiomers of catalyst needed for R/S. Bulky ketones slow. Enones can give 1,4.",
    mechanism="CBS oxazaborolidine coordinates BH₃ → B-H delivered to activated ketone → 6-membered TS → enantioselective hydride transfer.",
    refs="Corey, Bakshi, Shibata, JACS 1987, 109, 5551.",
    related="Noyori, Alpine-Borane, Midland, Luche",
    functional_groups="ketone, alcohol, chiral")

_nr("Shiina Macrolactonization",
    aliases=["Shiina macrolactonization", "Shiina esterification"],
    category="Functional Group Interconversion", type_="Macrolactonization",
    summary="MNBA (2-methyl-6-nitrobenzoic anhydride)-mediated macrolactonization, milder than Yamaguchi.",
    conditions="Seco-acid + MNBA, DMAP (cat.), Et₃N, CH₂Cl₂, high dilution (0.002 M), rt.",
    substrate="ω-Hydroxy acids → macrolactones (often with better yields than Yamaguchi for acid-sensitive substrates).",
    limitations="High dilution essential. MNBA must be freshly prepared.",
    mechanism="Seco-acid + MNBA → mixed anhydride → DMAP activation → intramolecular cyclization → macrolactone.",
    refs="Shiina, Kubota, Oshiumi, Hashizume, J. Org. Chem. 2004, 69, 1822.",
    related="Yamaguchi, Mitsunobu, Corey-Nicolaou",
    functional_groups="carboxylic acid, alcohol, macrolactone")

_nr("Suzuki-Miyaura Coupling",
    aliases=["Suzuki-Miyaura"],
    category="C-C Bond Formation", type_="Cross-Coupling",
    summary="Pd-catalyzed coupling of organoboron compounds with organic halides — the most widely used cross-coupling.",
    conditions="ArX + ArB(OH)₂, Pd(PPh₃)₄ or Pd(dppf)Cl₂, K₂CO₃ or K₃PO₄, dioxane/H₂O, 80 °C. Also: SPhos, XPhos ligands for challenging substrates.",
    substrate="Aryl/vinyl boronic acids/esters + aryl/vinyl halides/triflates → biaryls, dienes. Also alkyl-B (sp3-sp2).",
    limitations="Aryl chlorides need bulky ligands (SPhos). Boronic acid protodeborylation. Base required.",
    mechanism="Pd(0) + ArX → oxidative addition → Pd(II) + base → transmetalation with ArB(OH)₂ → reductive elimination → Ar-Ar + Pd(0).",
    refs="Miyaura, Suzuki, Chem. Rev. 1995, 95, 2457. Nobel Prize 2010.",
    related="Heck, Negishi, Stille, Kumada, Hiyama",
    functional_groups="aryl halide, boronic acid, biaryl")

_nr("Trost Asymmetric Allylic Alkylation",
    aliases=["Trost AAA", "asymmetric allylic alkylation"],
    category="C-C Bond Formation", type_="Asymmetric Allylation",
    summary="Pd-catalyzed AAA with Trost's DACH-phenyl ligand for enantioselective allylation.",
    conditions="Allyl acetate/carbonate + nucleophile (malonate, amine, phenol), Pd₂(dba)₃, Trost ligand, CH₂Cl₂, rt.",
    substrate="Allyl substrates + soft nucleophiles → branched/linear products with >90% ee.",
    limitations="Hard nucleophiles give poor ee. Linear/branched selectivity can be challenging.",
    mechanism="Pd(0) → oxidative addition → π-allyl-Pd(II) → chiral environment from Trost ligand → nucleophilic attack from differentiated face → product.",
    refs="Trost, Van Vranken, Chem. Rev. 1996, 96, 395.",
    related="Tsuji-Trost, Pd-allyl chemistry",
    functional_groups="allyl acetate, malonate, amine")

_nr("Swern Oxidation",
    aliases=["Swern oxidation"],
    category="Oxidation", type_="Alcohol Oxidation",
    summary="DMSO/oxalyl chloride-mediated oxidation of primary/secondary alcohols to aldehydes/ketones.",
    conditions="(COCl)₂, DMSO, CH₂Cl₂, –78 °C → add ROH → Et₃N → warm to rt.",
    substrate="1° alcohols → aldehydes (no over-oxidation). 2° alcohols → ketones. Very mild and selective.",
    limitations="Must maintain –78 °C during activation. Generates Me₂S (foul smell). Excess Et₃N needed.",
    mechanism="DMSO + (COCl)₂ → chlorosulfonium → ROH attacks S → alkoxysulfonium → Et₃N → intramolecular elimination → C=O.",
    refs="Omura, Swern, Tetrahedron 1978, 34, 1651.",
    related="Dess-Martin, PCC, TEMPO, Parikh-Doering",
    functional_groups="alcohol, aldehyde, ketone")

_nr("Barton-McCombie Deoxygenation",
    aliases=["Barton-McCombie", "radical deoxygenation"],
    category="Radical", type_="Deoxygenation",
    summary="Radical-mediated removal of OH groups: alcohol → xanthate/thiocarbamate → radical reduction → deoxygenated product.",
    conditions="ROH → xanthate (NaH, CS₂, MeI) → Bu₃SnH, AIBN, toluene, reflux. Modern: (TMS)₃SiH replaces tin.",
    substrate="Secondary and tertiary alcohols → deoxygenated. Tolerates many functional groups.",
    limitations="Tin waste (toxic). Xanthate preparation adds step. Primary alcohols less efficient.",
    mechanism="Bu₃Sn• abstracts S from xanthate → radical on carbon → β-scission → C radical → H from Bu₃SnH → product + new Bu₃Sn•.",
    refs="Barton, McCombie, J. Chem. Soc., Perkin Trans. 1 1975, 1574.",
    related="Barton decarboxylation",
    functional_groups="alcohol, xanthate, radical")

_nr("Corey-Winter Olefin Synthesis",
    aliases=["Corey-Winter"],
    category="C-C Bond Formation", type_="Olefination",
    summary="1,2-Diols → cyclic thionocarbonate → treatment with P(OMe)₃ → alkene (syn-elimination).",
    conditions="Diol + thiophosgene or thiocarbonyldiimidazole → cyclic thionocarbonate → P(OMe)₃, heat → alkene.",
    substrate="1,2-Diols → alkenes stereospecifically (syn-elimination: threo → Z, erythro → E).",
    limitations="Multi-step. Thiophosgene toxic. Better reagents (TCDI) available.",
    mechanism="P(OMe)₃ attacks S → betaine → cheletropic extrusion of CO/S → alkene (suprafacial, syn-periplanar).",
    refs="Corey, Winter, JACS 1963, 85, 2677.",
    related="Ramberg-Bäcklund, Bamford-Stevens",
    functional_groups="diol, alkene, thionocarbonate")

_nr("Mitsunobu Reaction",
    aliases=["Mitsunobu reaction"],
    category="Substitution", type_="Inversion SN2",
    summary="DIAD/PPh₃-mediated substitution of alcohols with inversion of stereochemistry.",
    conditions="ROH + NuH + DIAD (or DEAD) + PPh₃, THF, 0 °C → rt. Nu = carboxylates, phenols, phthalimides, sulfonamides, azides.",
    substrate="2° alcohols → substituted with inversion. Converts ROH → esters, ethers, azides, amines (Gabriel).",
    limitations="Atom-poor (2 eq reagents generate stoichiometric Ph₃P=O + DIAD-H₂). Purification of byproducts difficult.",
    mechanism="DIAD + PPh₃ → betaine → protonation by NuH → ion pair → SN2 on activated alcohol → product with inversion.",
    refs="Mitsunobu, Yamada, Bull. Chem. Soc. Jpn. 1967, 40, 2380.",
    related="Appel, Williamson, Gabriel",
    functional_groups="alcohol, ester, ether, inversion")

_nr("Grubbs Metathesis",
    aliases=["cross metathesis", "CM", "olefin cross metathesis"],
    category="C-C Bond Formation", type_="Metathesis",
    summary="Ru-catalyzed olefin cross-metathesis for convergent C=C bond formation between two alkenes.",
    conditions="Alkene A + alkene B, Grubbs II or Hoveyda-Grubbs, CH₂Cl₂, 40 °C. Ethylene removed (N₂ sparge).",
    substrate="Terminal alkenes + acrylates, styrenes, or allyl systems → cross products. Type I/II/III selectivity rules.",
    limitations="Homodimerization competes. Catalyst type-selectivity must be matched. Electron-rich olefins are Type I.",
    mechanism="Ru=CHR initiates → [2+2] cycloaddition → metallacyclobutane → retro-[2+2] → new Ru=CHR' → repeat → statistical/cross product.",
    refs="Chatterjee, Choi, Sanders, Grubbs, JACS 2003, 125, 11360.",
    related="RCM, ROMP, Olefin Metathesis",
    functional_groups="alkene, acrylate, styrene")

_nr("Oppolzer Sultam",
    aliases=["Oppolzer auxiliary", "camphorsultam auxiliary"],
    category="C-C Bond Formation", type_="Asymmetric Alkylation",
    summary="Bornane-10,2-sultam (Oppolzer's camphorsultam) as chiral auxiliary for asymmetric Diels-Alder, aldol, and alkylation.",
    conditions="N-Acylsultam + LDA → enolate → electrophile (alkyl halide, aldehyde, dienophile) → auxiliary removal (LiOH/H₂O₂).",
    substrate="Diels-Alder: dienophile on sultam → endo adduct with >98% de. Aldol: syn-selective.",
    limitations="Stoichiometric chiral auxiliary. Multi-step attachment/removal.",
    mechanism="Enolate chelated by SO₂ → electrophile approaches from less-hindered exo face → high de. Li or Ti enolates.",
    refs="Oppolzer, Tetrahedron 1987, 43, 1969.",
    related="Evans, Myers, Enders",
    functional_groups="acryloyl, aldehyde, Diels-Alder, sultam")

# =============================================================================
# Lookup Functions
# =============================================================================


def lookup_reaction(query: str) -> list[dict]:
    """Search named reactions by name, category, type, or substrate keyword."""
    q = query.lower().strip()
    results = []
    seen = set()
    # Exact match first
    if q in NAMED_REACTIONS:
        r = NAMED_REACTIONS[q]
        if r["name"] not in seen:
            results.append(r)
            seen.add(r["name"])

    # Fuzzy match: search across all fields
    for key, r in NAMED_REACTIONS.items():
        if r["name"] in seen:
            continue
        searchable = " ".join([
            r["name"], " ".join(r["aliases"]), r["category"], r["type"],
            r["summary"], r["conditions"], r["substrate"], r["functional_groups"],
        ]).lower()
        if q in searchable:
            results.append(r)
            seen.add(r["name"])

    return results[:15]


# =============================================================================
# Protecting Groups Database
# =============================================================================

# Stability scale: "S" = stable, "M" = moderately stable, "L" = labile
# Conditions: acid, base, nuc (nucleophile), ox (oxidation), red (reduction), h2_pd (hydrogenolysis)

PROTECTING_GROUPS: dict[str, list[dict]] = {
    "hydroxyl": [
        {"name": "TMS", "full_name": "Trimethylsilyl ether", "protection": "TMSCl, imidazole or Et3N, DMF or CH2Cl2", "deprotection": "TBAF, HF·py, AcOH/H2O/THF, or K2CO3/MeOH", "stability": {"acid": "L", "base": "M", "nuc": "M", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Least stable silyl ether. Removed under mildly acidic conditions."},
        {"name": "TBS (TBDMS)", "full_name": "tert-Butyldimethylsilyl ether", "protection": "TBSCl, imidazole, DMF (or TBSOTf, 2,6-lutidine, CH2Cl2)", "deprotection": "TBAF/THF, HF·py, AcOH/H2O, p-TsOH/MeOH", "stability": {"acid": "M", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Most popular silyl ether. ~10⁴× more stable than TMS. 1° OH selectivity over 2° with TBSCl/imid."},
        {"name": "TIPS", "full_name": "Triisopropylsilyl ether", "protection": "TIPSCl, imidazole, DMF or TIPSOTf, 2,6-lutidine", "deprotection": "TBAF, HF·py (slower than TBS)", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Very stable silyl ether. Steric bulk provides acid stability. ~10²× more stable than TBS toward acid."},
        {"name": "TES", "full_name": "Triethylsilyl ether", "protection": "TESCl, imidazole, DMF", "deprotection": "TBAF, AcOH, PPTS", "stability": {"acid": "L", "base": "S", "nuc": "M", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Intermediate stability. More stable than TMS, less than TBS."},
        {"name": "TBDPS", "full_name": "tert-Butyldiphenylsilyl ether", "protection": "TBDPSCl, imidazole, DMF", "deprotection": "TBAF (slow, sometimes requires heat), HF·py", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Very stable. UV-active (diphenyl). Selective for 1° OH. Compatible with many reaction conditions."},
        {"name": "Bn (Benzyl)", "full_name": "Benzyl ether", "protection": "BnBr, NaH, DMF or THF; or BnOC(=NH)CCl3, TfOH (Dudley)", "deprotection": "H2, Pd/C (hydrogenolysis); DDQ (PMB-selective); Li, NH3(l)", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "M", "red": "M", "h2_pd": "L"}, "notes": "Removed by hydrogenolysis. Orthogonal to silyl ethers. Very stable to acid/base."},
        {"name": "PMB", "full_name": "p-Methoxybenzyl ether", "protection": "PMBCl, NaH, DMF; or PMBTCA, BF3·Et2O", "deprotection": "DDQ, CH2Cl2/H2O; CAN; TFA; H2, Pd/C", "stability": {"acid": "M", "base": "S", "nuc": "S", "ox": "L", "red": "M", "h2_pd": "L"}, "notes": "Oxidatively labile (DDQ removes selectively in presence of Bn). More acid-labile than Bn."},
        {"name": "Ac (Acetyl)", "full_name": "Acetyl ester", "protection": "Ac2O, pyridine or DMAP; AcCl, Et3N", "deprotection": "K2CO3/MeOH (Zemplén); LiOH/THF-H2O; NaOMe/MeOH; NH3/MeOH", "stability": {"acid": "S", "base": "L", "nuc": "L", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Removed by mild base (transesterification in MeOH). Cheap, easy. Not stable to Grignard, LAH, or basic conditions."},
        {"name": "Piv (Pivaloyl)", "full_name": "Pivaloyl ester", "protection": "PivCl, Et3N, DMAP", "deprotection": "DIBAL-H, -78°C; LiAlH4; KOH(aq)", "stability": {"acid": "S", "base": "M", "nuc": "M", "ox": "S", "red": "M", "h2_pd": "S"}, "notes": "Sterically hindered ester. More resistant to base than Ac. Selective removal with DIBAL."},
        {"name": "MOM", "full_name": "Methoxymethyl ether", "protection": "MOMCl, iPr2NEt, CH2Cl2 (carcinogen caution); or dimethoxymethane, P2O5", "deprotection": "HCl/MeOH; conc. HCl/THF; PPTS/MeOH or iPrOH/heat; TMSBr", "stability": {"acid": "L", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Acid-labile acetal. MOMCl is a suspected carcinogen — handle with care. Stable to base, nucleophiles, organometallics."},
        {"name": "THP", "full_name": "Tetrahydropyranyl ether", "protection": "DHP (3,4-dihydro-2H-pyran), PPTS or TsOH, CH2Cl2", "deprotection": "PPTS/MeOH or EtOH; dilute HCl; AcOH/H2O/THF", "stability": {"acid": "L", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Acid-labile acetal. Creates new stereocenter (diastereomers). Cheap, easy to install. Can complicate NMR."},
        {"name": "Allyl", "full_name": "Allyl ether", "protection": "AllBr, NaH, DMF", "deprotection": "Pd(PPh3)4, morpholine or PhSiH3/Pd(PPh3)4; or Ir catalyst/H2 then acid", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "M", "red": "M", "h2_pd": "M"}, "notes": "Removed by Pd(0)-catalyzed allyl transfer. Orthogonal to acid/base-labile groups."},
    ],
    "amino": [
        {"name": "Boc", "full_name": "tert-Butyloxycarbonyl", "protection": "Boc2O, NaOH or Et3N or DMAP, THF/H2O", "deprotection": "TFA/CH2Cl2 (25-50%); 4M HCl in dioxane; TMSOTf/2,6-lutidine", "stability": {"acid": "L", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Acid-labile. Orthogonal to Fmoc (base-labile). Most common amino PG along with Fmoc and Cbz. Stable to hydrogenolysis, base, nucleophiles."},
        {"name": "Fmoc", "full_name": "9-Fluorenylmethyloxycarbonyl", "protection": "Fmoc-Cl or Fmoc-OSu, Na2CO3, dioxane/H2O", "deprotection": "Piperidine/DMF (20%); DBU; Et2NH", "stability": {"acid": "S", "base": "L", "nuc": "L", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Base-labile. Orthogonal to Boc. Standard in SPPS (solid-phase peptide synthesis). UV-active (monitoring deprotection)."},
        {"name": "Cbz (Z)", "full_name": "Benzyloxycarbonyl", "protection": "CbzCl (benzyl chloroformate), NaOH or NaHCO3, H2O/dioxane", "deprotection": "H2, Pd/C; HBr/AcOH; TMSI", "stability": {"acid": "M", "base": "S", "nuc": "S", "ox": "M", "red": "M", "h2_pd": "L"}, "notes": "Removed by hydrogenolysis. Orthogonal to Boc. Stable to mild acid and base."},
        {"name": "Alloc", "full_name": "Allyloxycarbonyl", "protection": "Alloc-Cl, pyridine", "deprotection": "Pd(PPh3)4, PhSiH3 or morpholine", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "S", "red": "M", "h2_pd": "M"}, "notes": "Removed by Pd(0). Orthogonal to both Boc and Fmoc. Used in complex peptide synthesis."},
        {"name": "Ns (Nosyl)", "full_name": "2-Nitrobenzenesulfonyl", "protection": "NsCl, Et3N, CH2Cl2", "deprotection": "PhSH, K2CO3, DMF; or mercaptoacetic acid", "stability": {"acid": "S", "base": "S", "nuc": "L", "ox": "S", "red": "M", "h2_pd": "M"}, "notes": "Removed by thiolysis (Fukuyama conditions). Activates nitrogen for alkylation (Mitsunobu on sulfonamide N)."},
        {"name": "Ts (Tosyl)", "full_name": "p-Toluenesulfonyl", "protection": "TsCl, pyridine or Et3N", "deprotection": "Na/naphthalenide; Mg/MeOH; SmI2; HBr/AcOH (harsh)", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "S", "red": "M", "h2_pd": "S"}, "notes": "Very stable — difficult to remove. Often considered semi-permanent. Electron-withdrawing on N."},
        {"name": "Troc", "full_name": "2,2,2-Trichloroethoxycarbonyl", "protection": "TrocCl, pyridine", "deprotection": "Zn, AcOH; Cd/Pb couple", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "S", "red": "L", "h2_pd": "S"}, "notes": "Removed reductively (Zn). Orthogonal to Boc, Fmoc, Cbz. Used in carbohydrate chemistry."},
        {"name": "Bn (N-Benzyl)", "full_name": "N-Benzyl", "protection": "BnBr, K2CO3, DMF; or reductive amination with PhCHO/NaBH4", "deprotection": "H2, Pd/C; Birch reduction; CAN (for PMP)", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "M", "red": "S", "h2_pd": "L"}, "notes": "Removed by hydrogenolysis. Simple and robust. N-PMB variant removed by CAN or DDQ."},
    ],
    "carbonyl": [
        {"name": "1,3-Dioxolane", "full_name": "Ethylene acetal/ketal", "protection": "HOCH2CH2OH, p-TsOH, toluene, Dean-Stark; or (MeO)2CH2, TMSOTf", "deprotection": "Aqueous acid (HCl, p-TsOH, Amberlyst); acetone/H+; I2/acetone", "stability": {"acid": "L", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Standard carbonyl protection. Stable to base, nucleophiles, metal hydrides (NaBH4, LiAlH4), Grignard, organolithium."},
        {"name": "1,3-Dioxane", "full_name": "Propylene acetal/ketal", "protection": "HOCH2CH2CH2OH, p-TsOH, toluene, Dean-Stark", "deprotection": "Aqueous acid", "stability": {"acid": "L", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Slightly more stable than 1,3-dioxolane. 6-membered ring."},
        {"name": "Dithiane", "full_name": "1,3-Dithiane (thioacetal)", "protection": "HSCH2CH2CH2SH, BF3·Et2O or I2; or HSCH2CH2SH for dithiolane", "deprotection": "HgCl2/CaCO3; NBS; MeI then hydrolysis; Selectfluor", "stability": {"acid": "S", "base": "S", "nuc": "S", "ox": "L", "red": "S", "h2_pd": "S"}, "notes": "Acid-stable unlike O-acetals. Also umpolung synthon (dithiane anion = acyl anion equivalent). Hg or oxidative deprotection."},
        {"name": "Dimethyl acetal", "full_name": "Dimethyl acetal", "protection": "MeOH, p-TsOH or trimethyl orthoformate, p-TsOH", "deprotection": "Aqueous acid, acetone/p-TsOH", "stability": {"acid": "L", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Open-chain acetal. Slightly more acid-labile than cyclic acetals."},
    ],
    "carboxyl": [
        {"name": "Methyl ester", "full_name": "Methyl ester", "protection": "CH2N2 (diazomethane); MeOH/H+; MeI/K2CO3", "deprotection": "LiOH/THF-H2O; NaOH; Me3SnOH; BBr3 (for phenyl methyl ethers too)", "stability": {"acid": "S", "base": "L", "nuc": "L", "ox": "S", "red": "M", "h2_pd": "S"}, "notes": "Simplest ester PG. Removed by saponification. Stable to acid. Not stable to strong nucleophiles or base."},
        {"name": "Ethyl ester", "full_name": "Ethyl ester", "protection": "EtOH/H+; EtI, K2CO3", "deprotection": "LiOH/THF-H2O; NaOH; KOH", "stability": {"acid": "S", "base": "L", "nuc": "L", "ox": "S", "red": "M", "h2_pd": "S"}, "notes": "Similar to methyl ester. Slightly slower saponification."},
        {"name": "tert-Butyl ester", "full_name": "tert-Butyl ester", "protection": "Isobutylene, H2SO4; Boc2O, DMAP (for amino acids); tBuOH, DCC", "deprotection": "TFA; 4M HCl/dioxane; TMSOTf/2,6-lutidine", "stability": {"acid": "L", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Acid-labile (like Boc). Stable to base and nucleophiles. Cleaved under same conditions as Boc."},
        {"name": "Benzyl ester", "full_name": "Benzyl ester", "protection": "BnBr, K2CO3, DMF; BnOH, DCC, DMAP; PhCH2OH, H+", "deprotection": "H2, Pd/C; HBr/AcOH; TMSI", "stability": {"acid": "M", "base": "M", "nuc": "M", "ox": "M", "red": "M", "h2_pd": "L"}, "notes": "Removed by hydrogenolysis. Orthogonal to tBu ester."},
        {"name": "Allyl ester", "full_name": "Allyl ester", "protection": "AllBr, K2CO3; allyl alcohol, DCC", "deprotection": "Pd(PPh3)4, morpholine or PhSiH3", "stability": {"acid": "S", "base": "M", "nuc": "M", "ox": "S", "red": "M", "h2_pd": "M"}, "notes": "Removed by Pd(0)-catalyzed allyl transfer. Orthogonal to Boc and Fmoc."},
        {"name": "TMSE", "full_name": "2-(Trimethylsilyl)ethyl ester", "protection": "TMSCH2CH2OH, DCC, DMAP", "deprotection": "TBAF, THF", "stability": {"acid": "M", "base": "S", "nuc": "S", "ox": "S", "red": "S", "h2_pd": "S"}, "notes": "Removed by fluoride (β-elimination of TMS). Very stable to base."},
    ],
}


def lookup_pg(query: str, functional_group: str | None = None) -> list[dict]:
    """Search protecting groups by name, functional group, or condition."""
    q = query.lower().strip()
    results = []

    groups_to_search = {}
    if functional_group:
        fg = functional_group.lower().strip()
        if fg in PROTECTING_GROUPS:
            groups_to_search[fg] = PROTECTING_GROUPS[fg]
    else:
        groups_to_search = PROTECTING_GROUPS

    for fg_name, pgs in groups_to_search.items():
        for pg in pgs:
            searchable = " ".join([
                pg["name"], pg.get("full_name", ""), pg.get("protection", ""),
                pg.get("deprotection", ""), pg.get("notes", ""),
            ]).lower()
            if q in searchable or q in pg["name"].lower():
                results.append({"functional_group": fg_name, **pg})

    return results


# =============================================================================
# Workup Procedures
# =============================================================================

WORKUP_PROCEDURES: dict[str, dict] = {
    "standard_aqueous": {
        "name": "Standard Aqueous Workup",
        "when": "Most reactions in organic solvents (CH2Cl2, EtOAc, Et2O, toluene)",
        "steps": [
            "Dilute reaction with organic solvent (EtOAc or CH2Cl2)",
            "Transfer to separatory funnel",
            "Wash with water (1×), then sat. NaHCO3 (1×) to remove acids",
            "Wash with brine (1×) to break emulsions and dry organic layer",
            "Dry organic layer over Na2SO4 or MgSO4 (15-30 min)",
            "Filter off drying agent",
            "Concentrate under reduced pressure (rotovap, water bath ≤40°C)",
        ],
        "tips": "Add brine if emulsion forms. If product is water-soluble, use back-extraction with organic solvent (3×). Check pH of aqueous layers.",
    },
    "acid_base_extraction": {
        "name": "Acid-Base Extraction",
        "when": "Separating acids, bases, and neutrals in a mixture",
        "steps": [
            "Dissolve mixture in Et2O or CH2Cl2",
            "Extract with 1M HCl (2×) — removes basic compounds (amines) into aqueous layer",
            "Extract organic layer with 1M NaOH (2×) — removes acidic compounds (carboxylic acids, phenols) into aqueous layer",
            "Remaining organic layer contains neutral compounds — dry over Na2SO4, filter, concentrate",
            "To recover bases: basify HCl extract with NaOH → extract with organic solvent",
            "To recover acids: acidify NaOH extract with HCl → extract with organic solvent",
        ],
        "tips": "For phenols vs carboxylic acids: extract first with NaHCO3 (removes only acids), then NaOH (removes phenols).",
    },
    "fieser_workup_lah": {
        "name": "Fieser Workup (LiAlH4 Reduction)",
        "when": "After LiAlH4 (LAH) reductions — quenching the aluminum salts",
        "steps": [
            "Cool reaction to 0°C",
            "For n grams of LAH used, add SLOWLY and DROPWISE:",
            "  n mL water (CAUTION: H2 evolution! Exothermic!)",
            "  n mL 15% NaOH(aq)",
            "  3n mL water",
            "Stir vigorously until white granular precipitate forms (30-60 min)",
            "Filter through Celite, wash filter cake with Et2O or EtOAc",
            "Dry filtrate over Na2SO4, filter, concentrate",
        ],
        "tips": "If precipitate is gelatinous instead of granular, add more NaOH. Alternatively, use Rochelle's salt (sat. sodium potassium tartrate) workup: quench into cold sat. Rochelle's salt, stir until layers separate. Na2SO4·10H2O method also works for small scale.",
    },
    "fieser_workup_dibal": {
        "name": "Workup for DIBAL-H Reductions",
        "when": "After DIBAL-H reductions (ester → aldehyde, etc.)",
        "steps": [
            "Cool to -78°C (or 0°C)",
            "Quench with sat. Rochelle's salt (sodium potassium tartrate) solution",
            "Warm to RT and stir vigorously for 1-2 hours until layers separate clearly",
            "Separate layers, extract aqueous with organic solvent (2-3×)",
            "Dry organic layer, filter, concentrate",
        ],
        "tips": "Alternative: quench with MeOH at -78°C, then add sat. Rochelle's salt. Celite filtration of aluminum salts if cloudy. Na/K tartrate chelates aluminum → clear phase separation.",
    },
    "oxidative_workup": {
        "name": "Peroxide/Oxidant Quench (After Hydroboration etc.)",
        "when": "After reactions using boranes, peroxides, or other oxidants",
        "steps": [
            "If using H2O2: ensure all peroxide is decomposed before workup",
            "Quench excess borane with careful addition of MeOH or water (gas evolution!)",
            "For hydroboration-oxidation: NaOH/H2O2 step then standard aqueous workup",
            "Test for peroxides with starch-iodide paper before concentrating",
            "Never concentrate to dryness when peroxides might be present!",
        ],
        "tips": "Na2SO3 or Na2S2O3 quenches excess peroxide/peracid. Always test before evaporation. Ethereal solvents can form peroxides.",
    },
    "metal_removal": {
        "name": "Transition Metal Removal",
        "when": "After Pd, Ru, Cu, or other transition metal-catalyzed reactions",
        "steps": [
            "For Pd: filter through Celite pad, wash with hot EtOAc",
            "For residual Pd: stir with SiliaMetS Thiol (silica-bound thiol scavenger) or activated charcoal",
            "For Ru (metathesis): DMSO treatment or PPh3 to decompose Ru, then filter through silica plug",
            "For Cu (CuAAC, Chan-Lam): wash with sat. NH4Cl or EDTA solution",
            "Column chromatography as final purification if metal contamination persists",
        ],
        "tips": "ICP-MS or ICP-OES for quantitative metal analysis in final product. Pd limit in pharmaceuticals: <10 ppm.",
    },
    "grignard_quench": {
        "name": "Grignard/Organolithium Quench",
        "when": "After Grignard or organolithium additions",
        "steps": [
            "Cool to 0°C",
            "Quench SLOWLY with sat. NH4Cl(aq) (drop by drop initially)",
            "Stir until clear bilayer forms (may need HCl if gel forms)",
            "Separate, extract aqueous layer with organic solvent (2-3×)",
            "Wash combined organics with brine, dry over Na2SO4, filter, concentrate",
        ],
        "tips": "For ketone products, sat. NH4Cl is mild enough to avoid epimerization. For acid-sensitive products, use sat. NaHCO3 instead. If thick emulsion: add more water, HCl until pH ~2, or filter through Celite.",
    },
    "swern_workup": {
        "name": "Swern Oxidation Workup",
        "when": "After Swern or related DMSO-based oxidations",
        "steps": [
            "Reaction mixture warms to RT after Et3N addition",
            "Dilute with CH2Cl2, wash with water (2×)",
            "Wash with 1M HCl (1×) to remove Et3N",
            "Wash with sat. NaHCO3 (1×), brine (1×)",
            "Dry over Na2SO4, filter, concentrate (fume hood — DMS smell!)",
        ],
        "tips": "Column chromatography usually needed. Me2S byproduct has strong odor — work in fume hood. Bleach trap for DMS waste.",
    },
}


def lookup_workup(query: str) -> list[dict]:
    """Search workup procedures by name or keyword."""
    q = query.lower().strip()
    results = []
    for key, proc in WORKUP_PROCEDURES.items():
        searchable = " ".join([proc["name"], proc["when"], key, " ".join(proc["steps"]), proc.get("tips", "")]).lower()
        if q in searchable or q in key:
            results.append(proc)
    return results


# =============================================================================
# Solvent Guide
# =============================================================================

SOLVENTS: list[dict] = [
    # Common organic solvents: name, bp, density, polarity_index, dielectric_constant, miscible_with_water, common_uses
    {"name": "Pentane", "bp": 36, "density": 0.626, "polarity_index": 0.0, "dielectric": 1.84, "water_misc": False, "uses": "Column chromatography (non-polar), extractions", "safety": "Highly flammable, low flash point"},
    {"name": "Hexane(s)", "bp": 69, "density": 0.659, "polarity_index": 0.1, "dielectric": 1.88, "water_misc": False, "uses": "Column chromatography, recrystallization solvent pair with EtOAc", "safety": "Flammable, neurotoxic (n-hexane)"},
    {"name": "Heptane", "bp": 98, "density": 0.684, "polarity_index": 0.1, "dielectric": 1.92, "water_misc": False, "uses": "Column chromatography alternative to hexane, crystallization", "safety": "Flammable"},
    {"name": "Cyclohexane", "bp": 81, "density": 0.779, "polarity_index": 0.2, "dielectric": 2.02, "water_misc": False, "uses": "Crystallization, NMR solvent (C6D12)", "safety": "Flammable"},
    {"name": "Petroleum ether", "bp": "40-60", "density": 0.64, "polarity_index": 0.1, "dielectric": 1.9, "water_misc": False, "uses": "Column chromatography, extraction", "safety": "Highly flammable, mixture of hydrocarbons"},
    {"name": "Diethyl ether", "bp": 35, "density": 0.713, "polarity_index": 2.8, "dielectric": 4.34, "water_misc": False, "uses": "Grignard reactions, extractions, crystallization. Forms peroxides!", "safety": "Extremely flammable, peroxide-forming, low bp"},
    {"name": "MTBE (methyl tert-butyl ether)", "bp": 55, "density": 0.740, "polarity_index": 2.5, "dielectric": 2.6, "water_misc": False, "uses": "Extraction solvent (safer than Et2O), less prone to peroxides", "safety": "Flammable, but safer than Et2O"},
    {"name": "THF", "bp": 66, "density": 0.889, "polarity_index": 4.0, "dielectric": 7.58, "water_misc": True, "uses": "Organometallic reactions, Grignard, LDA, universal solvent. Inhibited with BHT.", "safety": "Flammable, peroxide-forming. Distill from Na/benzophenone for dry THF."},
    {"name": "2-MeTHF", "bp": 80, "density": 0.854, "polarity_index": 3.5, "dielectric": 6.97, "water_misc": False, "uses": "Green alternative to THF. Better phase separation. Organometallics.", "safety": "Flammable. Bio-derived."},
    {"name": "1,4-Dioxane", "bp": 101, "density": 1.034, "polarity_index": 4.8, "dielectric": 2.25, "water_misc": True, "uses": "Pd cross-coupling, Buchwald-Hartwig. Good for high-temp reactions.", "safety": "Peroxide-forming, suspected carcinogen"},
    {"name": "Dichloromethane (DCM)", "bp": 40, "density": 1.325, "polarity_index": 3.1, "dielectric": 8.93, "water_misc": False, "uses": "Extractions (heavier than water!), column chromatography, reactions at RT", "safety": "Not flammable but toxic (possible carcinogen). Denser than water."},
    {"name": "Chloroform", "bp": 61, "density": 1.483, "polarity_index": 4.1, "dielectric": 4.81, "water_misc": False, "uses": "NMR solvent (CDCl3), extractions", "safety": "Toxic, possible carcinogen. Contains EtOH stabilizer. Denser than water."},
    {"name": "1,2-Dichloroethane (DCE)", "bp": 83, "density": 1.253, "polarity_index": 3.5, "dielectric": 10.36, "water_misc": False, "uses": "Lewis acid reactions (higher bp alternative to DCM)", "safety": "Toxic, carcinogenic"},
    {"name": "Ethyl acetate (EtOAc)", "bp": 77, "density": 0.902, "polarity_index": 4.4, "dielectric": 6.02, "water_misc": False, "uses": "Column chromatography, extraction, recrystallization", "safety": "Flammable, pleasant smell"},
    {"name": "Acetone", "bp": 56, "density": 0.791, "polarity_index": 5.1, "dielectric": 20.7, "water_misc": True, "uses": "Cleaning glassware, Jones oxidation, dry ice baths, fast-evaporating", "safety": "Highly flammable"},
    {"name": "Acetonitrile (MeCN)", "bp": 82, "density": 0.786, "polarity_index": 5.8, "dielectric": 37.5, "water_misc": True, "uses": "HPLC solvent, reactions with electrophilic reagents, Ritter reaction", "safety": "Flammable, toxic"},
    {"name": "Methanol (MeOH)", "bp": 65, "density": 0.791, "polarity_index": 5.1, "dielectric": 32.7, "water_misc": True, "uses": "Recrystallization, hydrogenolysis, Zemplén deacetylation", "safety": "Flammable, toxic (blindness)"},
    {"name": "Ethanol (EtOH)", "bp": 78, "density": 0.789, "polarity_index": 5.2, "dielectric": 24.5, "water_misc": True, "uses": "Recrystallization, green solvent, reactions", "safety": "Flammable"},
    {"name": "Isopropanol (iPrOH)", "bp": 82, "density": 0.786, "polarity_index": 3.9, "dielectric": 17.9, "water_misc": True, "uses": "Cleaning, MPV reduction, crystallization", "safety": "Flammable"},
    {"name": "DMF", "bp": 153, "density": 0.944, "polarity_index": 6.4, "dielectric": 36.7, "water_misc": True, "uses": "SN2 reactions, cross-coupling, peptide synthesis. Dissolves many salts.", "safety": "Reproductive toxin. Difficult to remove. Decomposes to Me2NH + CO at high temp with base."},
    {"name": "DMA (N,N-dimethylacetamide)", "bp": 165, "density": 0.937, "polarity_index": 6.5, "dielectric": 37.8, "water_misc": True, "uses": "Alternative to DMF. Similar properties but slightly higher bp.", "safety": "Reproductive toxin"},
    {"name": "NMP", "bp": 204, "density": 1.028, "polarity_index": 6.7, "dielectric": 32.2, "water_misc": True, "uses": "High-temp reactions, Stille coupling, polymer dissolution", "safety": "Reproductive toxin"},
    {"name": "DMSO", "bp": 189, "density": 1.100, "polarity_index": 7.2, "dielectric": 46.7, "water_misc": True, "uses": "Swern oxidation, Kornblum oxidation, high-temp SNAr, DMSO-d6 NMR", "safety": "Penetrates skin rapidly (carries dissolved chemicals through skin!)"},
    {"name": "Toluene", "bp": 111, "density": 0.867, "polarity_index": 2.4, "dielectric": 2.38, "water_misc": False, "uses": "Azeotropic drying (Dean-Stark), Grubbs metathesis, general non-polar solvent", "safety": "Flammable, toxic"},
    {"name": "Benzene", "bp": 80, "density": 0.879, "polarity_index": 2.7, "dielectric": 2.28, "water_misc": False, "uses": "Legacy solvent (replaced by toluene). NMR (C6D6)", "safety": "Known carcinogen — avoid use!"},
    {"name": "Xylene(s)", "bp": "138-144", "density": 0.860, "polarity_index": 2.5, "dielectric": 2.37, "water_misc": False, "uses": "High-temp reactions, histology", "safety": "Flammable, toxic"},
    {"name": "Water", "bp": 100, "density": 1.000, "polarity_index": 10.2, "dielectric": 80.1, "water_misc": True, "uses": "Green solvent, aqueous reactions, on-water chemistry", "safety": "Non-toxic. Incompatible with many reagents."},
    {"name": "Pyridine", "bp": 115, "density": 0.982, "polarity_index": 5.3, "dielectric": 12.4, "water_misc": True, "uses": "Base/solvent for acylations (Ac2O/py), nucleophilic catalysis, chromium oxidations", "safety": "Toxic, unpleasant smell, flammable"},
    {"name": "Triethylamine (Et3N)", "bp": 89, "density": 0.726, "polarity_index": 1.8, "dielectric": 2.42, "water_misc": False, "uses": "Base (non-nucleophilic), Hünig's base alternative, scavenger for HX", "safety": "Flammable, corrosive, strong smell"},
    {"name": "DIPEA (Hünig's base)", "bp": 127, "density": 0.742, "polarity_index": 2.0, "dielectric": 3.4, "water_misc": False, "uses": "Non-nucleophilic base, peptide coupling, cross-coupling", "safety": "Flammable"},
    {"name": "HFIP (hexafluoroisopropanol)", "bp": 59, "density": 1.596, "polarity_index": 5.3, "dielectric": 16.7, "water_misc": True, "uses": "Strong H-bond donor, metal-free oxidations, hypervalent iodine chemistry", "safety": "Expensive, toxic"},
    {"name": "TFE (trifluoroethanol)", "bp": 74, "density": 1.383, "polarity_index": 4.3, "dielectric": 26.7, "water_misc": True, "uses": "Acid-surrogate solvent, C-H activation, radical reactions", "safety": "Toxic"},
]


# Chromatography solvent systems
CHROM_SOLVENT_SYSTEMS = [
    {"system": "Hexanes / EtOAc", "range": "Non-polar to moderate", "notes": "Most common system. Start at 5% EtOAc, increase. Standard for general organic compounds."},
    {"system": "Hexanes / Et2O", "range": "Non-polar to moderate", "notes": "Weaker eluent than EtOAc at same %. Good for separating hydrocarbons and halides."},
    {"system": "Hexanes / Acetone", "range": "Non-polar to moderate", "notes": "Stronger eluent than EtOAc at same %. Good for polar compounds."},
    {"system": "CH2Cl2 / MeOH", "range": "Moderate to very polar", "notes": "For polar compounds (amines, acids, amides). Start at 1-2% MeOH. Add 0.1% Et3N for basic compounds or 0.1% AcOH for acids."},
    {"system": "CH2Cl2 / EtOAc", "range": "Moderate polarity", "notes": "Good intermediate polarity range. Useful when hex/EtOAc gradient is too shallow."},
    {"system": "Toluene / EtOAc", "range": "Non-polar to moderate", "notes": "Alternative to hex/EtOAc. Better for UV-active compounds."},
    {"system": "EtOAc / MeOH", "range": "Polar to very polar", "notes": "Very polar compounds. Use sparingly — can dissolve silica."},
    {"system": "CHCl3 / MeOH / NH4OH", "range": "Very polar (amines)", "notes": "For very polar amines. Ratio: 90:9:1 to 80:18:2."},
]


def lookup_solvent(query: str) -> list[dict]:
    """Search solvents by name, property, or use."""
    q = query.lower().strip()
    results = []
    for s in SOLVENTS:
        searchable = " ".join([s["name"], s["uses"], s.get("safety", "")]).lower()
        if q in searchable:
            results.append(s)
    return results


# =============================================================================
# Cooling Baths
# =============================================================================

COOLING_BATHS: list[dict] = [
    {"temp": 0, "recipe": "Ice / water", "notes": "Most common. Maintain ice:water ~2:1. Replace ice as needed."},
    {"temp": -5, "recipe": "Ice / NaCl (rock salt)", "notes": "~1:3 salt:ice by mass."},
    {"temp": -10, "recipe": "Ice / NaCl (rock salt)", "notes": "~1:3 salt:ice. Pack well."},
    {"temp": -15, "recipe": "Ice / CaCl2 (anhydrous, 1:2.5 by mass)", "notes": "Calcium chloride/ice gives lower temp than NaCl/ice."},
    {"temp": -20, "recipe": "Ice / NaCl (1:3 by mass)", "notes": "Or commercial -20°C freezer."},
    {"temp": -40, "recipe": "Dry ice / acetonitrile", "notes": "Acetonitrile/dry ice slurry. bp(MeCN) = 82°C. Maintains -40°C."},
    {"temp": -42, "recipe": "Dry ice / acetonitrile", "notes": "Acetonitrile/CO2 slurry."},
    {"temp": -45, "recipe": "Dry ice / acetonitrile or dry ice / CH3CN", "notes": "Practical range -40 to -45°C."},
    {"temp": -78, "recipe": "Dry ice / acetone", "notes": "THE standard low-temp bath. Sublimation temp of CO2. Add dry ice chunks to acetone. Replenish dry ice as needed."},
    {"temp": -78, "recipe": "Dry ice / isopropanol", "notes": "Better thermal contact than dry ice/acetone. Same temp."},
    {"temp": -84, "recipe": "Ethyl acetate / liquid N2", "notes": "EtOAc slush bath. Add LN2 slowly to EtOAc."},
    {"temp": -94, "recipe": "Toluene / liquid N2", "notes": "Toluene slush. Add LN2 carefully to toluene."},
    {"temp": -98, "recipe": "Methanol / liquid N2", "notes": "MeOH slush. Used for extremely low temp operations."},
    {"temp": -116, "recipe": "Ethanol / liquid N2", "notes": "EtOH slush. Very low temperature."},
    {"temp": -131, "recipe": "Pentane / liquid N2", "notes": "Pentane slush. Near LN2 temperatures."},
    {"temp": -196, "recipe": "Liquid nitrogen", "notes": "Boiling point of N2. Use Dewar flask. Extremely cold — frostbite danger!"},
    # Heat sources
    {"temp": 40, "recipe": "Water bath", "notes": "Simple and precise. Use with rotovap."},
    {"temp": 60, "recipe": "Water bath", "notes": "Good for gentle heating."},
    {"temp": 80, "recipe": "Water bath or oil bath", "notes": "Near boiling water. Transition to oil bath."},
    {"temp": 100, "recipe": "Boiling water bath or oil bath", "notes": "Steam bath for gentle reflux of many solvents."},
    {"temp": 120, "recipe": "Oil bath (silicone or mineral oil)", "notes": "Standard oil bath range. Monitor with thermometer."},
    {"temp": 150, "recipe": "Oil bath (silicone oil)", "notes": "Mid-range oil bath. Silicone oil to ~250°C."},
    {"temp": 200, "recipe": "Oil bath or sand bath", "notes": "High-temp oil bath or sand bath. Mineral oil limit ~250°C."},
    {"temp": 250, "recipe": "Sand bath or aluminum block", "notes": "Sand or metal heating block for high temp."},
]


def lookup_cooling_bath(target_temp: float | None = None) -> list[dict]:
    """Find cooling/heating bath for target temperature. Returns closest options."""
    if target_temp is None:
        return COOLING_BATHS

    # Sort by distance from target
    sorted_baths = sorted(COOLING_BATHS, key=lambda b: abs(b["temp"] - target_temp))
    return sorted_baths[:5]


# =============================================================================
# TLC Stains
# =============================================================================

TLC_STAINS: list[dict] = [
    {
        "name": "UV (254 nm)",
        "detects": "Aromatic compounds, conjugated systems, any UV-active chromophore",
        "recipe": "Use UV-active TLC plates (F254). Shine 254 nm UV lamp — spots appear dark on green background.",
        "notes": "First method to try. Non-destructive. Does not detect saturated aliphatics. Short wavelength (254 nm) quenches plate fluorescence. Long wavelength (365 nm) for fluorescent compounds."
    },
    {
        "name": "KMnO4 (Potassium permanganate)",
        "detects": "Alkenes, alkynes, aldehydes, 1°/2° alcohols, most oxidizable compounds. UNIVERSAL stain.",
        "recipe": "1.5 g KMnO4, 10 g K2CO3, 1.25 mL 10% NaOH, 200 mL water. Dip TLC plate, heat with heat gun.",
        "notes": "Yellow spots on purple background. Most universal stain — use as default if unsure. Does NOT stain tertiary alcohols, alkyl halides, or non-oxidizable compounds."
    },
    {
        "name": "p-Anisaldehyde",
        "detects": "Almost everything: alcohols, terpenes, steroids, sugars, nucleosides. Multi-colored spots!",
        "recipe": "135 mL EtOH + 5 mL conc. H2SO4 + 1.5 mL AcOH + 3.7 mL p-anisaldehyde. Dip, heat strongly.",
        "notes": "Multi-colored stain (blue, purple, green, pink, orange) — very useful for distinguishing compounds. Keep refrigerated. Particularly good for terpenes, carbohydrates."
    },
    {
        "name": "CAM (Ceric ammonium molybdate / Hanessian's stain)",
        "detects": "Universal: alcohols, amines, ethers, most organics. Especially good for alcohol/amine mixtures.",
        "recipe": "5 g ammonium molybdate + 2 g ceric sulfate in 200 mL 10% H2SO4. Dip, heat.",
        "notes": "Dark blue/black spots on light blue background. Versatile and sensitive. Colors develop on heating."
    },
    {
        "name": "PMA (Phosphomolybdic acid)",
        "detects": "Universal: almost anything reducing. Especially lipids, steroids.",
        "recipe": "10% phosphomolybdic acid in EtOH. Dip, heat.",
        "notes": "Dark green/blue spots on yellow-green background. Very sensitive. Works for almost everything. Charring at high temp."
    },
    {
        "name": "Vanillin",
        "detects": "Alcohols, aldehydes, ketones, terpenes, terpenoids, higher MW compounds. Multi-colored.",
        "recipe": "15 g vanillin + 250 mL EtOH + 2.5 mL conc. H2SO4. Dip, heat.",
        "notes": "Multi-colored spots. Similar utility to anisaldehyde. Good for terpenes and steroids."
    },
    {
        "name": "Ninhydrin",
        "detects": "Primary amines, amino acids. Purple/pink spots ('Ruhemann's purple').",
        "recipe": "0.3 g ninhydrin in 100 mL n-butanol + 3 mL AcOH. Dip, heat at 110°C.",
        "notes": "Specific for amines. Proline/hydroxyproline give yellow. Secondary amines give brown/orange. Standard amino acid detection."
    },
    {
        "name": "Dragendorff's reagent",
        "detects": "Alkaloids, heterocyclic nitrogen compounds, quaternary ammonium compounds.",
        "recipe": "Solution A: 0.85 g bismuth subnitrate in 10 mL AcOH + 40 mL H2O. Solution B: 8 g KI in 20 mL H2O. Mix A + B, add 20 mL AcOH + 100 mL H2O.",
        "notes": "Orange/red spots on yellow background. Specific for nitrogen heterocycles."
    },
    {
        "name": "DNP (2,4-Dinitrophenylhydrazine)",
        "detects": "Aldehydes and ketones specifically.",
        "recipe": "0.4 g 2,4-DNP in 2 mL H2SO4 + 30 mL EtOH + 10 mL H2O. Dip (no heating needed).",
        "notes": "Orange/red spots. Highly specific for carbonyls. Classic qualitative test. Also used for derivatization for melting point ID."
    },
    {
        "name": "Bromocresol green",
        "detects": "Carboxylic acids and other acidic compounds.",
        "recipe": "0.04 g bromocresol green in 100 mL EtOH + enough NaOH to make blue. Dip (no heating).",
        "notes": "Yellow spots on blue background. Specific for acids. Non-destructive (at RT)."
    },
    {
        "name": "Iodine chamber",
        "detects": "Unsaturated compounds, aromatics, many organics. Reversible stain.",
        "recipe": "Place a few crystals of I2 in a sealed chamber (jar with lid). Place TLC plate inside, wait.",
        "notes": "Brown spots. Reversible — spots fade after removal from chamber. Non-destructive. Good for quick check."
    },
    {
        "name": "Seebach's stain ('Magic stain')",
        "detects": "Universal stain. Detects virtually all organic compounds.",
        "recipe": "2.5 g phosphomolybdic acid + 1 g ceric sulfate + 6 mL conc. H2SO4 + 94 mL H2O. Dip, heat.",
        "notes": "Multi-colored spots on green background. Combines PMA + CAN. Extremely versatile."
    },
    {
        "name": "FeCl3 (Ferric chloride)",
        "detects": "Phenols, enols, hydroxamic acids, sulfhydryl groups.",
        "recipe": "1% FeCl3 in water or EtOH. Dip or spray.",
        "notes": "Purple/blue/green/red spots depending on phenol structure. Classic phenol test."
    },
]


def lookup_tlc_stain(query: str) -> list[dict]:
    """Search TLC stains by name, functional group detected, or keyword."""
    q = query.lower().strip()
    results = []
    for stain in TLC_STAINS:
        searchable = " ".join([stain["name"], stain["detects"], stain.get("notes", "")]).lower()
        if q in searchable:
            results.append(stain)
    # If no match, return all as reference
    if not results and q in ("all", "list", "help", ""):
        return TLC_STAINS
    return results


# =============================================================================
# Column Chromatography Guide
# =============================================================================

COLUMN_GUIDE = {
    "still_guidelines": {
        "title": "Still's Flash Chromatography Guidelines (Still, Kahn, Mitra, JOC 1978)",
        "rules": [
            "Target Rf of 0.2-0.35 in your chosen solvent system (by TLC)",
            "If Rf < 0.2: increase polarity; if Rf > 0.35: decrease polarity",
            "Column diameter: use ~1 inch per 200 mg crude for easy separations",
            "Silica height: 6 inches (15 cm) of silica for standard separation",
            "Loading: ≤100:1 silica/crude (by mass) for difficult separations (ΔRf < 0.1)",
            "Loading: ~30:1 silica/crude for easy separations (ΔRf > 0.2)",
            "Apply sample as concentrated solution or dry-loaded on silica/Celite",
            "Elute with ~10 column volumes of solvent",
            "Collect fractions: ~5 mL per cm of column diameter",
            "Use air pressure (flash) to maintain flow rate ~2 inches/min",
        ],
    },
    "troubleshooting": [
        {"problem": "Bands are streaking/tailing", "causes": "Overloaded column, compound too polar for system, acidic/basic compounds adsorbing to silica", "solutions": "Use less crude, increase polarity, add 0.1-1% Et3N (for amines) or 0.1-1% AcOH (for acids) to mobile phase. Try Florisil or alumina instead of silica."},
        {"problem": "No separation (bands overlap)", "causes": "Rf values too similar, inadequate column length, wrong solvent system", "solutions": "Try different solvent systems (e.g., switch from hex/EtOAc to hex/acetone). Use longer column (more silica). Use smaller fractions."},
        {"problem": "Product stuck on column", "causes": "Compound too polar, strong adsorption", "solutions": "Flush with high-polarity solvent (MeOH, then EtOAc+1%Et3N). Use reverse-phase or alumina."},
        {"problem": "Fronting (leading edge)", "causes": "Sample overload, solvent mismatch", "solutions": "Load less material. Ensure loading solvent matches or is less polar than mobile phase."},
        {"problem": "Poor recovery / decomposition on column", "causes": "Acid-sensitive compounds on silica (pH ~5-6), long contact time", "solutions": "Use deactivated silica (1-5% Et3N pretreated), basic alumina, or Florisil. Run column faster. Use MPLC."},
        {"problem": "Baseline impurity co-eluting", "causes": "Silica dissolution, BHT from solvents, plasticizer from tubing", "solutions": "Use distilled solvents. Avoid PVC tubing. Pre-rinse column with mobile phase."},
    ],
    "column_sizes": [
        {"diameter_cm": 1.0, "crude_mg": "40-100", "silica_g": "4-10", "fraction_ml": "3-5"},
        {"diameter_cm": 1.5, "crude_mg": "100-200", "silica_g": "10-20", "fraction_ml": "5-8"},
        {"diameter_cm": 2.0, "crude_mg": "200-400", "silica_g": "20-40", "fraction_ml": "8-12"},
        {"diameter_cm": 2.5, "crude_mg": "400-800", "silica_g": "40-80", "fraction_ml": "12-15"},
        {"diameter_cm": 3.0, "crude_mg": "800-1200", "silica_g": "60-120", "fraction_ml": "15-20"},
        {"diameter_cm": 4.0, "crude_mg": "1200-2000", "silica_g": "100-200", "fraction_ml": "20-30"},
        {"diameter_cm": 5.0, "crude_mg": "2000-5000", "silica_g": "200-500", "fraction_ml": "30-50"},
    ],
    "rf_rules": [
        "Rf = distance traveled by compound / distance traveled by solvent front",
        "Rf is reproducible only if: same silica, same chamber saturation, same temperature",
        "Rf 0.0-0.1: too polar for this system — increase mobile phase polarity",
        "Rf 0.1-0.2: acceptable but slow elution — consider more polar system",
        "Rf 0.2-0.35: IDEAL range for column chromatography",
        "Rf 0.35-0.5: OK for easy separations or gradient",
        "Rf 0.5-1.0: compound too non-polar for this system — decrease polarity",
        "If two spots have ΔRf ≥ 0.15: easy separation (~20:1 silica:crude)",
        "If ΔRf ≈ 0.1: moderate separation (~50:1 silica:crude)",
        "If ΔRf < 0.05: difficult separation (consider different system, gradient, or prep HPLC)",
    ],
    "dry_loading": {
        "description": "Dry loading (preferred for best resolution)",
        "steps": [
            "Dissolve crude in minimum volume of CH2Cl2 or EtOAc",
            "Add ~2× mass Celite or silica gel (relative to crude)",
            "Evaporate solvent on rotovap → free-flowing powder",
            "Add powder on top of packed column (below sand layer)",
            "This gives a narrow initial band → better resolution",
        ],
    },
}


def lookup_column_guide(query: str) -> dict | list:
    """Search column chromatography guide by topic."""
    q = query.lower().strip()
    if "still" in q or "guideline" in q or "rule" in q or "general" in q:
        return COLUMN_GUIDE["still_guidelines"]
    if "troubleshoot" in q or "problem" in q or "streak" in q or "tail" in q:
        return COLUMN_GUIDE["troubleshooting"]
    if "size" in q or "diameter" in q or "how much silica" in q:
        return COLUMN_GUIDE["column_sizes"]
    if "rf" in q or "retention" in q:
        return COLUMN_GUIDE["rf_rules"]
    if "dry load" in q or "loading" in q:
        return COLUMN_GUIDE["dry_loading"]
    if "solvent" in q or "system" in q or "eluent" in q:
        return CHROM_SOLVENT_SYSTEMS
    # Return everything
    return COLUMN_GUIDE


# =============================================================================
# Buffer Recipes — 25 common laboratory buffers
# =============================================================================

BUFFER_RECIPES: list[dict] = [
    {
        "name": "PBS (Phosphate-Buffered Saline)",
        "pH": "7.4",
        "recipe": "137 mM NaCl, 2.7 mM KCl, 10 mM Na₂HPO₄, 1.8 mM KH₂PO₄",
        "preparation": "Dissolve 8.0 g NaCl, 0.2 g KCl, 1.44 g Na₂HPO₄, 0.24 g KH₂PO₄ in 800 mL dH₂O. Adjust pH to 7.4 with HCl. Add water to 1 L. Autoclave.",
        "category": "physiological",
        "pKa": 7.20,
        "range": "6.2–8.2",
        "notes": "Most common biological buffer. Do NOT use with Ca²⁺/Mg²⁺-sensitive assays (phosphate chelates divalent cations).",
    },
    {
        "name": "10× PBS",
        "pH": "7.4",
        "recipe": "1.37 M NaCl, 27 mM KCl, 100 mM Na₂HPO₄, 18 mM KH₂PO₄",
        "preparation": "Dissolve 80 g NaCl, 2 g KCl, 14.4 g Na₂HPO₄, 2.4 g KH₂PO₄ in 800 mL dH₂O. Adjust pH. Add water to 1 L.",
        "category": "physiological",
        "pKa": 7.20,
        "range": "6.2–8.2",
        "notes": "Dilute 1:10 before use. Check pH after dilution.",
    },
    {
        "name": "TBS (Tris-Buffered Saline)",
        "pH": "7.6",
        "recipe": "20 mM Tris-HCl, 150 mM NaCl",
        "preparation": "Dissolve 2.42 g Tris base + 8.77 g NaCl in 900 mL dH₂O. Adjust pH to 7.6 with conc. HCl (~3.2 mL). Add water to 1 L.",
        "category": "physiological",
        "pKa": 8.07,
        "range": "7.0–9.0",
        "notes": "Alternative to PBS when phosphate interferes. Tris pKa is temperature-sensitive: ΔpKa/°C = −0.031.",
    },
    {
        "name": "Tris-HCl",
        "pH": "7.4–8.8",
        "recipe": "50 mM Tris-HCl",
        "preparation": "Dissolve 6.06 g Tris base in 800 mL dH₂O. Adjust pH with conc. HCl at the working temperature. Add water to 1 L.",
        "category": "general",
        "pKa": 8.07,
        "range": "7.0–9.0",
        "notes": "Adjust pH at working temperature (pKa drops ~0.03/°C). Common for protein purification, electrophoresis.",
    },
    {
        "name": "HEPES",
        "pH": "7.2–7.6",
        "recipe": "25 mM HEPES",
        "preparation": "Dissolve 5.96 g HEPES in 900 mL dH₂O. Adjust pH with NaOH. Add water to 1 L.",
        "category": "biological",
        "pKa": 7.55,
        "range": "6.8–8.2",
        "notes": "Good's buffer. Minimal metal binding. Preferred for cell culture. Does not interfere with BCA assay.",
    },
    {
        "name": "MES",
        "pH": "5.5–6.7",
        "recipe": "50 mM MES",
        "preparation": "Dissolve 9.76 g MES in 900 mL dH₂O. Adjust pH with NaOH. Add water to 1 L.",
        "category": "biological",
        "pKa": 6.15,
        "range": "5.5–6.7",
        "notes": "Good's buffer for slightly acidic conditions. Useful for ion-exchange chromatography.",
    },
    {
        "name": "MOPS",
        "pH": "6.5–7.9",
        "recipe": "20 mM MOPS",
        "preparation": "Dissolve 4.19 g MOPS in 900 mL dH₂O. Adjust pH with NaOH. Add water to 1 L.",
        "category": "biological",
        "pKa": 7.20,
        "range": "6.5–7.9",
        "notes": "Good's buffer. Common for RNA gel electrophoresis (MOPS buffer = 20 mM MOPS, 5 mM NaOAc, 1 mM EDTA, pH 7.0).",
    },
    {
        "name": "PIPES",
        "pH": "6.1–7.5",
        "recipe": "25 mM PIPES",
        "preparation": "Dissolve 7.56 g PIPES in 900 mL dH₂O. May need slight warming to dissolve. Adjust pH with NaOH. Add water to 1 L.",
        "category": "biological",
        "pKa": 6.76,
        "range": "6.1–7.5",
        "notes": "Good's buffer. Low metal binding. Relatively insoluble in water — dissolve with NaOH.",
    },
    {
        "name": "Citrate buffer",
        "pH": "3.0–6.2",
        "recipe": "0.1 M citric acid / 0.1 M sodium citrate",
        "preparation": "Mix 0.1 M citric acid (21.01 g/L) with 0.1 M trisodium citrate dihydrate (29.41 g/L) to desired pH.",
        "category": "general",
        "pKa": 3.13,
        "range": "3.0–6.2",
        "notes": "Triprotic: pKa₁=3.13, pKa₂=4.76, pKa₃=6.40. Chelates metals. Common for protein crystallization.",
    },
    {
        "name": "Acetate buffer",
        "pH": "3.6–5.6",
        "recipe": "0.1 M acetic acid / 0.1 M sodium acetate",
        "preparation": "Mix 0.1 M acetic acid (5.72 mL glacial AcOH/L) with 0.1 M NaOAc (8.20 g/L) to desired pH.",
        "category": "general",
        "pKa": 4.76,
        "range": "3.6–5.6",
        "notes": "Volatile — good for mass spectrometry. Common for protein purification at acidic pH.",
    },
    {
        "name": "Phosphate buffer (Sørensen)",
        "pH": "5.8–8.0",
        "recipe": "Varies — mix monobasic + dibasic sodium phosphate",
        "preparation": "Mix 0.2 M NaH₂PO₄ (27.60 g/L) with 0.2 M Na₂HPO₄ (28.39 g/L) in ratios from tables to desired pH.",
        "category": "general",
        "pKa": 7.20,
        "range": "5.8–8.0",
        "notes": "Autoclavable. Avoid with proteins sensitive to phosphate. Precipitates with Ca²⁺, Mg²⁺, Mn²⁺, Fe³⁺.",
    },
    {
        "name": "Carbonate-bicarbonate buffer",
        "pH": "9.2–10.8",
        "recipe": "0.1 M Na₂CO₃ / 0.1 M NaHCO₃",
        "preparation": "Mix 0.1 M Na₂CO₃ (10.60 g/L) with 0.1 M NaHCO₃ (8.40 g/L) to desired pH.",
        "category": "general",
        "pKa": 10.33,
        "range": "9.2–10.8",
        "notes": "High pH buffer. Common for ELISA coating. CO₂ loss shifts pH — keep sealed.",
    },
    {
        "name": "Glycine-HCl",
        "pH": "2.2–3.6",
        "recipe": "0.1 M glycine adjusted with HCl",
        "preparation": "Dissolve 7.51 g glycine in 900 mL dH₂O. Adjust pH with conc. HCl. Add water to 1 L.",
        "category": "general",
        "pKa": 2.35,
        "range": "2.2–3.6",
        "notes": "Common for protein elution from affinity columns (low pH elution).",
    },
    {
        "name": "TAE (Tris-Acetate-EDTA)",
        "pH": "8.0",
        "recipe": "40 mM Tris, 20 mM acetic acid, 1 mM EDTA",
        "preparation": "50×: 242 g Tris base + 57.1 mL glacial acetic acid + 100 mL 0.5 M EDTA pH 8.0 → 1 L. Dilute 1:50.",
        "category": "electrophoresis",
        "pKa": 8.07,
        "range": "7.5–8.5",
        "notes": "DNA gel electrophoresis. Better resolution than TBE but lower buffering capacity.",
    },
    {
        "name": "TBE (Tris-Borate-EDTA)",
        "pH": "8.3",
        "recipe": "89 mM Tris, 89 mM boric acid, 2 mM EDTA",
        "preparation": "10×: 108 g Tris base + 55 g boric acid + 40 mL 0.5 M EDTA pH 8.0 → 1 L. Dilute 1:10.",
        "category": "electrophoresis",
        "pKa": 8.07,
        "range": "8.0–8.5",
        "notes": "DNA/RNA gel electrophoresis. Higher buffering capacity than TAE. Better for small fragments (<1 kb).",
    },
    {
        "name": "SDS-PAGE running buffer (Laemmli)",
        "pH": "8.3",
        "recipe": "25 mM Tris, 192 mM glycine, 0.1% SDS",
        "preparation": "10×: 30.3 g Tris base + 144 g glycine + 10 g SDS → 1 L dH₂O. Dilute 1:10. Do NOT adjust pH.",
        "category": "electrophoresis",
        "pKa": 8.07,
        "range": "8.3",
        "notes": "Tris-glycine system for SDS-PAGE. pH is set by the buffer combination — do NOT add acid/base.",
    },
    {
        "name": "Western transfer buffer (Towbin)",
        "pH": "8.3",
        "recipe": "25 mM Tris, 192 mM glycine, 20% methanol",
        "preparation": "3.03 g Tris base + 14.4 g glycine + 200 mL methanol → 1 L dH₂O. Do NOT adjust pH.",
        "category": "electrophoresis",
        "pKa": 8.07,
        "range": "8.3",
        "notes": "For Western blot transfer. Methanol helps protein binding to membrane. Reduce to 10% for large proteins (>100 kDa).",
    },
    {
        "name": "TE buffer",
        "pH": "8.0",
        "recipe": "10 mM Tris-HCl pH 8.0, 1 mM EDTA",
        "preparation": "Mix 1 mL 1 M Tris-HCl pH 8.0 + 0.2 mL 0.5 M EDTA pH 8.0 → 100 mL dH₂O. Autoclave.",
        "category": "nucleic acid",
        "pKa": 8.07,
        "range": "7.5–8.5",
        "notes": "Standard DNA/RNA storage buffer. EDTA chelates nucleases. Use TE-low (0.1 mM EDTA) for PCR templates.",
    },
    {
        "name": "RIPA lysis buffer",
        "pH": "7.4",
        "recipe": "50 mM Tris-HCl pH 7.4, 150 mM NaCl, 1% NP-40, 0.5% sodium deoxycholate, 0.1% SDS",
        "preparation": "Mix components in dH₂O. Add protease inhibitors fresh before use.",
        "category": "lysis",
        "pKa": 8.07,
        "range": "7.0–8.0",
        "notes": "Strong lysis for total protein extraction. Add protease + phosphatase inhibitors fresh. Keep on ice.",
    },
    {
        "name": "Lysis buffer (mild, NP-40)",
        "pH": "7.4",
        "recipe": "50 mM Tris-HCl pH 7.4, 150 mM NaCl, 1% NP-40",
        "preparation": "Mix components in dH₂O. Add protease inhibitors fresh before use.",
        "category": "lysis",
        "pKa": 8.07,
        "range": "7.0–8.0",
        "notes": "Mild lysis preserving protein-protein interactions. Good for co-IP.",
    },
]


def lookup_buffer(query: str) -> list[dict]:
    """Search buffer recipes by name, pH range, category, or keyword."""
    q = query.lower().strip()
    results = []
    for buf in BUFFER_RECIPES:
        searchable = f"{buf['name']} {buf.get('category', '')} {buf.get('notes', '')} {buf.get('recipe', '')}".lower()
        if q in searchable:
            results.append(buf)
    # If no exact match, try individual words
    if not results:
        words = q.split()
        for buf in BUFFER_RECIPES:
            searchable = f"{buf['name']} {buf.get('category', '')} {buf.get('notes', '')}".lower()
            if all(w in searchable for w in words):
                results.append(buf)
    # Try pH range
    if not results:
        try:
            target_ph = float(q)
            for buf in BUFFER_RECIPES:
                try:
                    low, high = buf["range"].split("–")
                    if float(low) <= target_ph <= float(high):
                        results.append(buf)
                except (ValueError, KeyError):
                    pass
        except ValueError:
            pass
    return results


# =============================================================================
# Amino Acid Properties — 20 canonical + common modifications
# =============================================================================

AMINO_ACIDS: list[dict] = [
    {"code1": "G", "code3": "Gly", "name": "Glycine",       "mw": 57.02,  "pKa_carboxyl": 2.35, "pKa_amino": 9.78, "pKa_side": None,  "pI": 5.97, "hydropathy": -0.4,  "class": "nonpolar",   "notes": "Smallest AA. Achiral. Flexible — disrupts helices."},
    {"code1": "A", "code3": "Ala", "name": "Alanine",       "mw": 71.04,  "pKa_carboxyl": 2.34, "pKa_amino": 9.87, "pKa_side": None,  "pI": 6.01, "hydropathy": 1.8,   "class": "nonpolar",   "notes": "Strong helix former. Reference for hydrophobicity scales."},
    {"code1": "V", "code3": "Val", "name": "Valine",        "mw": 99.07,  "pKa_carboxyl": 2.29, "pKa_amino": 9.74, "pKa_side": None,  "pI": 5.97, "hydropathy": 4.2,   "class": "nonpolar",   "notes": "Branched-chain AA. Beta-branched — disfavors helix."},
    {"code1": "L", "code3": "Leu", "name": "Leucine",       "mw": 113.08, "pKa_carboxyl": 2.33, "pKa_amino": 9.74, "pKa_side": None,  "pI": 5.98, "hydropathy": 3.8,   "class": "nonpolar",   "notes": "Most abundant AA in proteins. Strong helix former."},
    {"code1": "I", "code3": "Ile", "name": "Isoleucine",    "mw": 113.08, "pKa_carboxyl": 2.32, "pKa_amino": 9.76, "pKa_side": None,  "pI": 6.02, "hydropathy": 4.5,   "class": "nonpolar",   "notes": "Beta-branched. Two chiral centers. Most hydrophobic canonical AA."},
    {"code1": "P", "code3": "Pro", "name": "Proline",       "mw": 97.05,  "pKa_carboxyl": 1.95, "pKa_amino": 10.64, "pKa_side": None, "pI": 6.48, "hydropathy": -1.6,  "class": "nonpolar",   "notes": "Cyclic (imino acid). Helix breaker. cis-trans isomerism. Rigid."},
    {"code1": "F", "code3": "Phe", "name": "Phenylalanine", "mw": 147.07, "pKa_carboxyl": 2.20, "pKa_amino": 9.31, "pKa_side": None,  "pI": 5.49, "hydropathy": 2.8,   "class": "aromatic",   "notes": "UV absorbs weakly at 257 nm. Aromatic π-stacking."},
    {"code1": "W", "code3": "Trp", "name": "Tryptophan",    "mw": 186.08, "pKa_carboxyl": 2.46, "pKa_amino": 9.41, "pKa_side": None,  "pI": 5.89, "hydropathy": -0.9,  "class": "aromatic",   "notes": "Largest AA. UV absorbs at 280 nm (ε=5500). Rare in proteins. Fluorescent (λex=295, λem=348 nm)."},
    {"code1": "M", "code3": "Met", "name": "Methionine",    "mw": 131.04, "pKa_carboxyl": 2.13, "pKa_amino": 9.28, "pKa_side": None,  "pI": 5.74, "hydropathy": 1.9,   "class": "nonpolar",   "notes": "Start codon (AUG). Thioether side chain. Easily oxidized to sulfoxide."},
    {"code1": "S", "code3": "Ser", "name": "Serine",        "mw": 87.03,  "pKa_carboxyl": 2.19, "pKa_amino": 9.21, "pKa_side": None,  "pI": 5.68, "hydropathy": -0.8,  "class": "polar uncharged", "notes": "Phosphorylation site. Active site nucleophile (serine proteases)."},
    {"code1": "T", "code3": "Thr", "name": "Threonine",     "mw": 101.05, "pKa_carboxyl": 2.09, "pKa_amino": 9.10, "pKa_side": None,  "pI": 5.60, "hydropathy": -0.7,  "class": "polar uncharged", "notes": "Phosphorylation site. Beta-branched. Two chiral centers."},
    {"code1": "C", "code3": "Cys", "name": "Cysteine",      "mw": 103.01, "pKa_carboxyl": 1.92, "pKa_amino": 10.70, "pKa_side": 8.37, "pI": 5.07, "hydropathy": 2.5,   "class": "polar uncharged", "notes": "Disulfide bonds (S-S). Strong nucleophile (thiol). UV: ε_280=125 (disulfide). Sensitive to oxidation."},
    {"code1": "Y", "code3": "Tyr", "name": "Tyrosine",      "mw": 163.06, "pKa_carboxyl": 2.20, "pKa_amino": 9.21, "pKa_side": 10.46, "pI": 5.66, "hydropathy": -1.3, "class": "aromatic",   "notes": "Phosphorylation site. UV absorbs at 280 nm (ε=1490). Phenol OH."},
    {"code1": "H", "code3": "His", "name": "Histidine",     "mw": 137.06, "pKa_carboxyl": 1.78, "pKa_amino": 9.18, "pKa_side": 6.04,  "pI": 7.59, "hydropathy": -3.2,  "class": "positively charged", "notes": "Imidazole (pKa ≈ 6). Protonated at physiological pH ~10%. Acid-base catalysis. Metal coordination (Zn, Cu, Fe)."},
    {"code1": "K", "code3": "Lys", "name": "Lysine",        "mw": 128.09, "pKa_carboxyl": 2.16, "pKa_amino": 9.06, "pKa_side": 10.54, "pI": 9.74, "hydropathy": -3.9,  "class": "positively charged", "notes": "Charged at physiological pH. Ubiquitination, acetylation, methylation target."},
    {"code1": "R", "code3": "Arg", "name": "Arginine",      "mw": 156.10, "pKa_carboxyl": 1.82, "pKa_amino": 8.99, "pKa_side": 12.48, "pI": 10.76, "hydropathy": -4.5, "class": "positively charged", "notes": "Guanidinium — always charged at physiological pH. Strongest base among AAs. Salt bridges."},
    {"code1": "D", "code3": "Asp", "name": "Aspartate",     "mw": 115.03, "pKa_carboxyl": 1.99, "pKa_amino": 9.90, "pKa_side": 3.90,  "pI": 2.77, "hydropathy": -3.5,  "class": "negatively charged", "notes": "Carboxylate side chain. Metal coordination. Succinimide formation with adjacent Asn/Ser."},
    {"code1": "E", "code3": "Glu", "name": "Glutamate",     "mw": 129.04, "pKa_carboxyl": 2.10, "pKa_amino": 9.47, "pKa_side": 4.07,  "pI": 3.22, "hydropathy": -3.5,  "class": "negatively charged", "notes": "Carboxylate. Neurotransmitter. Strong helix former. Longer side chain than Asp."},
    {"code1": "N", "code3": "Asn", "name": "Asparagine",    "mw": 114.04, "pKa_carboxyl": 2.14, "pKa_amino": 8.72, "pKa_side": None,  "pI": 5.41, "hydropathy": -3.5,  "class": "polar uncharged", "notes": "Amide side chain. N-glycosylation site (Asn-X-Ser/Thr). Deamidation-prone."},
    {"code1": "Q", "code3": "Gln", "name": "Glutamine",     "mw": 128.06, "pKa_carboxyl": 2.17, "pKa_amino": 9.13, "pKa_side": None,  "pI": 5.65, "hydropathy": -3.5,  "class": "polar uncharged", "notes": "Amide side chain. Deamidation-prone (→ Glu). Polyglutamine expansions in disease."},
]

# Build lookup dicts
_AA_BY_CODE1 = {aa["code1"]: aa for aa in AMINO_ACIDS}
_AA_BY_CODE3 = {aa["code3"].lower(): aa for aa in AMINO_ACIDS}
_AA_BY_NAME = {aa["name"].lower(): aa for aa in AMINO_ACIDS}


def lookup_amino_acid(query: str) -> list[dict]:
    """Look up amino acid by 1-letter code, 3-letter code, name, or property class."""
    q = query.strip()

    # Single character → 1-letter code
    if len(q) == 1:
        aa = _AA_BY_CODE1.get(q.upper())
        return [aa] if aa else []

    # 3-letter code
    aa = _AA_BY_CODE3.get(q.lower())
    if aa:
        return [aa]

    # Full name
    aa = _AA_BY_NAME.get(q.lower())
    if aa:
        return [aa]

    # Class or property search
    q_lower = q.lower()
    results = []
    for aa in AMINO_ACIDS:
        searchable = f"{aa['name']} {aa['class']} {aa.get('notes', '')}".lower()
        if q_lower in searchable:
            results.append(aa)

    # If "all" or empty — return everything
    if not results and q_lower in ("all", "table", "list", ""):
        return AMINO_ACIDS

    return results


# =============================================================================
# NMR Solvent Reference — residual solvent peaks for common NMR solvents
# =============================================================================

NMR_SOLVENTS: list[dict] = [
    {"name": "Chloroform-d", "formula": "CDCl₃", "abbrev": "CDCl3",
     "h_residual": 7.26, "c_residual": 77.16, "c_multiplicity": "triplet",
     "bp": 61, "density": 1.492, "water_peak": 1.56,
     "notes": "Most common NMR solvent. Dissolves most organic compounds. Slightly acidic — avoid with base-sensitive compounds."},
    {"name": "DMSO-d₆", "formula": "(CD₃)₂SO", "abbrev": "DMSO-d6",
     "h_residual": 2.50, "c_residual": 39.52, "c_multiplicity": "septet",
     "bp": 189, "density": 1.190, "water_peak": 3.33,
     "notes": "Universal solvent. Hygroscopic. High bp — hard to remove. Dissolves salts, peptides, polymers."},
    {"name": "Methanol-d₄", "formula": "CD₃OD", "abbrev": "MeOD",
     "h_residual": 3.31, "c_residual": 49.00, "c_multiplicity": "septet",
     "bp": 65, "density": 0.888, "water_peak": 4.87,
     "notes": "Good for polar compounds. Exchangeable protons (OH, NH) disappear. Two residual peaks: 3.31 (CHD₂) and 4.87 (OH)."},
    {"name": "D₂O", "formula": "D₂O", "abbrev": "D2O",
     "h_residual": 4.79, "c_residual": None, "c_multiplicity": None,
     "bp": 101, "density": 1.107, "water_peak": 4.79,
     "notes": "For water-soluble compounds. No ¹³C reference — use external standard (DSS, TSP). Exchangeable protons disappear. pH-dependent HDO peak."},
    {"name": "Acetone-d₆", "formula": "(CD₃)₂CO", "abbrev": "acetone-d6",
     "h_residual": 2.05, "c_residual": 29.84, "c_multiplicity": "septet",
     "bp": 56, "density": 0.872, "water_peak": 2.84,
     "notes": "Good for moderately polar compounds. Carbonyl ¹³C at 206.26 ppm. Low viscosity → sharp peaks."},
    {"name": "Acetonitrile-d₃", "formula": "CD₃CN", "abbrev": "MeCN-d3",
     "h_residual": 1.94, "c_residual": 1.32, "c_multiplicity": "septet",
     "bp": 82, "density": 0.844, "water_peak": 2.13,
     "notes": "Good for small molecules. CN carbon at 118.26 ppm. Low viscosity."},
    {"name": "Benzene-d₆", "formula": "C₆D₆", "abbrev": "C6D6",
     "h_residual": 7.16, "c_residual": 128.06, "c_multiplicity": "triplet",
     "bp": 80, "density": 0.950, "water_peak": 0.40,
     "notes": "Aromatic solvent. Shifts aromatic peaks upfield (ASIS effect). Carcinogen — handle with care."},
    {"name": "THF-d₈", "formula": "C₄D₈O", "abbrev": "THF-d8",
     "h_residual": "1.72, 3.58", "c_residual": "25.31, 67.21", "c_multiplicity": "quintets",
     "bp": 66, "density": 0.985, "water_peak": 2.46,
     "notes": "Good for organometallics, air-sensitive chemistry. Two sets of peaks (α and β CH₂)."},
    {"name": "DCM-d₂", "formula": "CD₂Cl₂", "abbrev": "CD2Cl2",
     "h_residual": 5.32, "c_residual": 53.84, "c_multiplicity": "quintet",
     "bp": 40, "density": 1.350, "water_peak": 1.52,
     "notes": "Similar to CDCl₃ but lower bp. Good for low-temperature NMR. Less acidic than CDCl₃."},
    {"name": "DMF-d₇", "formula": "(CD₃)₂NCDO", "abbrev": "DMF-d7",
     "h_residual": "2.75, 2.92, 8.03", "c_residual": "29.76, 34.89, 162.62", "c_multiplicity": "various",
     "bp": 153, "density": 1.044, "water_peak": 3.49,
     "notes": "Polar aprotic. Three residual ¹H peaks. Good for peptides, salts, organometallics."},
    {"name": "Pyridine-d₅", "formula": "C₅D₅N", "abbrev": "py-d5",
     "h_residual": "7.22, 7.58, 8.74", "c_residual": "123.87, 135.91, 149.90", "c_multiplicity": "triplets",
     "bp": 115, "density": 1.046, "water_peak": 4.95,
     "notes": "Basic solvent. Dissolves carbohydrates, nucleosides. Three ¹H residual peaks."},
    {"name": "TFA-d", "formula": "CF₃COOD", "abbrev": "TFA-d",
     "h_residual": 11.50, "c_residual": "116.6, 164.2", "c_multiplicity": "quartets (¹⁹F coupling)",
     "bp": 72, "density": 1.485, "water_peak": None,
     "notes": "Strong acid solvent. ¹³C shows quartets from ¹⁹F coupling. Dissolves peptides, proteins."},
]


def lookup_nmr_solvent(query: str) -> list[dict]:
    """Search NMR solvent reference by name or abbreviation."""
    q = query.lower().strip()
    results = []
    for sol in NMR_SOLVENTS:
        searchable = f"{sol['name']} {sol['abbrev']} {sol['formula']} {sol.get('notes', '')}".lower()
        if q in searchable:
            results.append(sol)
    if not results and q in ("all", "table", "list", ""):
        return NMR_SOLVENTS
    return results
