"""
labmate-mcp chemistry utilities module.

Pure-computation tools:
  - Isotope pattern calculation from molecular formula
  - CAS registry number validation
  - Scientific unit conversion
  - Periodic table element lookup
  - Buffer pH / Henderson-Hasselbalch calculations

No external API calls — all offline, zero latency.
"""

from __future__ import annotations

import math
import re
from itertools import product as iter_product
from typing import Any

# =============================================================================
# Periodic table — element data (Z ≤ 118)
# =============================================================================

ELEMENTS: dict[str, dict[str, Any]] = {
    "H":  {"Z": 1,  "name": "Hydrogen",     "mass": 1.008,    "group": 1,  "period": 1, "block": "s",  "en": 2.20, "config": "1s1",           "category": "nonmetal"},
    "He": {"Z": 2,  "name": "Helium",       "mass": 4.003,    "group": 18, "period": 1, "block": "s",  "en": None, "config": "1s2",           "category": "noble gas"},
    "Li": {"Z": 3,  "name": "Lithium",      "mass": 6.941,    "group": 1,  "period": 2, "block": "s",  "en": 0.98, "config": "[He] 2s1",      "category": "alkali metal"},
    "Be": {"Z": 4,  "name": "Beryllium",    "mass": 9.012,    "group": 2,  "period": 2, "block": "s",  "en": 1.57, "config": "[He] 2s2",      "category": "alkaline earth"},
    "B":  {"Z": 5,  "name": "Boron",        "mass": 10.81,    "group": 13, "period": 2, "block": "p",  "en": 2.04, "config": "[He] 2s2 2p1",  "category": "metalloid"},
    "C":  {"Z": 6,  "name": "Carbon",       "mass": 12.011,   "group": 14, "period": 2, "block": "p",  "en": 2.55, "config": "[He] 2s2 2p2",  "category": "nonmetal"},
    "N":  {"Z": 7,  "name": "Nitrogen",     "mass": 14.007,   "group": 15, "period": 2, "block": "p",  "en": 3.04, "config": "[He] 2s2 2p3",  "category": "nonmetal"},
    "O":  {"Z": 8,  "name": "Oxygen",       "mass": 15.999,   "group": 16, "period": 2, "block": "p",  "en": 3.44, "config": "[He] 2s2 2p4",  "category": "nonmetal"},
    "F":  {"Z": 9,  "name": "Fluorine",     "mass": 18.998,   "group": 17, "period": 2, "block": "p",  "en": 3.98, "config": "[He] 2s2 2p5",  "category": "halogen"},
    "Ne": {"Z": 10, "name": "Neon",         "mass": 20.180,   "group": 18, "period": 2, "block": "p",  "en": None, "config": "[He] 2s2 2p6",  "category": "noble gas"},
    "Na": {"Z": 11, "name": "Sodium",       "mass": 22.990,   "group": 1,  "period": 3, "block": "s",  "en": 0.93, "config": "[Ne] 3s1",      "category": "alkali metal"},
    "Mg": {"Z": 12, "name": "Magnesium",    "mass": 24.305,   "group": 2,  "period": 3, "block": "s",  "en": 1.31, "config": "[Ne] 3s2",      "category": "alkaline earth"},
    "Al": {"Z": 13, "name": "Aluminium",    "mass": 26.982,   "group": 13, "period": 3, "block": "p",  "en": 1.61, "config": "[Ne] 3s2 3p1",  "category": "post-transition metal"},
    "Si": {"Z": 14, "name": "Silicon",      "mass": 28.086,   "group": 14, "period": 3, "block": "p",  "en": 1.90, "config": "[Ne] 3s2 3p2",  "category": "metalloid"},
    "P":  {"Z": 15, "name": "Phosphorus",   "mass": 30.974,   "group": 15, "period": 3, "block": "p",  "en": 2.19, "config": "[Ne] 3s2 3p3",  "category": "nonmetal"},
    "S":  {"Z": 16, "name": "Sulfur",       "mass": 32.06,    "group": 16, "period": 3, "block": "p",  "en": 2.58, "config": "[Ne] 3s2 3p4",  "category": "nonmetal"},
    "Cl": {"Z": 17, "name": "Chlorine",     "mass": 35.45,    "group": 17, "period": 3, "block": "p",  "en": 3.16, "config": "[Ne] 3s2 3p5",  "category": "halogen"},
    "Ar": {"Z": 18, "name": "Argon",        "mass": 39.948,   "group": 18, "period": 3, "block": "p",  "en": None, "config": "[Ne] 3s2 3p6",  "category": "noble gas"},
    "K":  {"Z": 19, "name": "Potassium",    "mass": 39.098,   "group": 1,  "period": 4, "block": "s",  "en": 0.82, "config": "[Ar] 4s1",      "category": "alkali metal"},
    "Ca": {"Z": 20, "name": "Calcium",      "mass": 40.078,   "group": 2,  "period": 4, "block": "s",  "en": 1.00, "config": "[Ar] 4s2",      "category": "alkaline earth"},
    "Sc": {"Z": 21, "name": "Scandium",     "mass": 44.956,   "group": 3,  "period": 4, "block": "d",  "en": 1.36, "config": "[Ar] 3d1 4s2",  "category": "transition metal"},
    "Ti": {"Z": 22, "name": "Titanium",     "mass": 47.867,   "group": 4,  "period": 4, "block": "d",  "en": 1.54, "config": "[Ar] 3d2 4s2",  "category": "transition metal"},
    "V":  {"Z": 23, "name": "Vanadium",     "mass": 50.942,   "group": 5,  "period": 4, "block": "d",  "en": 1.63, "config": "[Ar] 3d3 4s2",  "category": "transition metal"},
    "Cr": {"Z": 24, "name": "Chromium",     "mass": 51.996,   "group": 6,  "period": 4, "block": "d",  "en": 1.66, "config": "[Ar] 3d5 4s1",  "category": "transition metal"},
    "Mn": {"Z": 25, "name": "Manganese",    "mass": 54.938,   "group": 7,  "period": 4, "block": "d",  "en": 1.55, "config": "[Ar] 3d5 4s2",  "category": "transition metal"},
    "Fe": {"Z": 26, "name": "Iron",         "mass": 55.845,   "group": 8,  "period": 4, "block": "d",  "en": 1.83, "config": "[Ar] 3d6 4s2",  "category": "transition metal"},
    "Co": {"Z": 27, "name": "Cobalt",       "mass": 58.933,   "group": 9,  "period": 4, "block": "d",  "en": 1.88, "config": "[Ar] 3d7 4s2",  "category": "transition metal"},
    "Ni": {"Z": 28, "name": "Nickel",       "mass": 58.693,   "group": 10, "period": 4, "block": "d",  "en": 1.91, "config": "[Ar] 3d8 4s2",  "category": "transition metal"},
    "Cu": {"Z": 29, "name": "Copper",       "mass": 63.546,   "group": 11, "period": 4, "block": "d",  "en": 1.90, "config": "[Ar] 3d10 4s1", "category": "transition metal"},
    "Zn": {"Z": 30, "name": "Zinc",         "mass": 65.38,    "group": 12, "period": 4, "block": "d",  "en": 1.65, "config": "[Ar] 3d10 4s2", "category": "transition metal"},
    "Ga": {"Z": 31, "name": "Gallium",      "mass": 69.723,   "group": 13, "period": 4, "block": "p",  "en": 1.81, "config": "[Ar] 3d10 4s2 4p1", "category": "post-transition metal"},
    "Ge": {"Z": 32, "name": "Germanium",    "mass": 72.63,    "group": 14, "period": 4, "block": "p",  "en": 2.01, "config": "[Ar] 3d10 4s2 4p2", "category": "metalloid"},
    "As": {"Z": 33, "name": "Arsenic",      "mass": 74.922,   "group": 15, "period": 4, "block": "p",  "en": 2.18, "config": "[Ar] 3d10 4s2 4p3", "category": "metalloid"},
    "Se": {"Z": 34, "name": "Selenium",     "mass": 78.96,    "group": 16, "period": 4, "block": "p",  "en": 2.55, "config": "[Ar] 3d10 4s2 4p4", "category": "nonmetal"},
    "Br": {"Z": 35, "name": "Bromine",      "mass": 79.904,   "group": 17, "period": 4, "block": "p",  "en": 2.96, "config": "[Ar] 3d10 4s2 4p5", "category": "halogen"},
    "Kr": {"Z": 36, "name": "Krypton",      "mass": 83.798,   "group": 18, "period": 4, "block": "p",  "en": 3.00, "config": "[Ar] 3d10 4s2 4p6", "category": "noble gas"},
    "Rb": {"Z": 37, "name": "Rubidium",     "mass": 85.468,   "group": 1,  "period": 5, "block": "s",  "en": 0.82, "config": "[Kr] 5s1",      "category": "alkali metal"},
    "Sr": {"Z": 38, "name": "Strontium",    "mass": 87.62,    "group": 2,  "period": 5, "block": "s",  "en": 0.95, "config": "[Kr] 5s2",      "category": "alkaline earth"},
    "Y":  {"Z": 39, "name": "Yttrium",      "mass": 88.906,   "group": 3,  "period": 5, "block": "d",  "en": 1.22, "config": "[Kr] 4d1 5s2",  "category": "transition metal"},
    "Zr": {"Z": 40, "name": "Zirconium",    "mass": 91.224,   "group": 4,  "period": 5, "block": "d",  "en": 1.33, "config": "[Kr] 4d2 5s2",  "category": "transition metal"},
    "Nb": {"Z": 41, "name": "Niobium",      "mass": 92.906,   "group": 5,  "period": 5, "block": "d",  "en": 1.60, "config": "[Kr] 4d4 5s1",  "category": "transition metal"},
    "Mo": {"Z": 42, "name": "Molybdenum",   "mass": 95.96,    "group": 6,  "period": 5, "block": "d",  "en": 2.16, "config": "[Kr] 4d5 5s1",  "category": "transition metal"},
    "Ru": {"Z": 44, "name": "Ruthenium",    "mass": 101.07,   "group": 8,  "period": 5, "block": "d",  "en": 2.20, "config": "[Kr] 4d7 5s1",  "category": "transition metal"},
    "Rh": {"Z": 45, "name": "Rhodium",      "mass": 102.906,  "group": 9,  "period": 5, "block": "d",  "en": 2.28, "config": "[Kr] 4d8 5s1",  "category": "transition metal"},
    "Pd": {"Z": 46, "name": "Palladium",    "mass": 106.42,   "group": 10, "period": 5, "block": "d",  "en": 2.20, "config": "[Kr] 4d10",     "category": "transition metal"},
    "Ag": {"Z": 47, "name": "Silver",       "mass": 107.868,  "group": 11, "period": 5, "block": "d",  "en": 1.93, "config": "[Kr] 4d10 5s1", "category": "transition metal"},
    "Cd": {"Z": 48, "name": "Cadmium",      "mass": 112.411,  "group": 12, "period": 5, "block": "d",  "en": 1.69, "config": "[Kr] 4d10 5s2", "category": "transition metal"},
    "In": {"Z": 49, "name": "Indium",       "mass": 114.818,  "group": 13, "period": 5, "block": "p",  "en": 1.78, "config": "[Kr] 4d10 5s2 5p1", "category": "post-transition metal"},
    "Sn": {"Z": 50, "name": "Tin",          "mass": 118.710,  "group": 14, "period": 5, "block": "p",  "en": 1.96, "config": "[Kr] 4d10 5s2 5p2", "category": "post-transition metal"},
    "Sb": {"Z": 51, "name": "Antimony",     "mass": 121.760,  "group": 15, "period": 5, "block": "p",  "en": 2.05, "config": "[Kr] 4d10 5s2 5p3", "category": "metalloid"},
    "Te": {"Z": 52, "name": "Tellurium",    "mass": 127.60,   "group": 16, "period": 5, "block": "p",  "en": 2.10, "config": "[Kr] 4d10 5s2 5p4", "category": "metalloid"},
    "I":  {"Z": 53, "name": "Iodine",       "mass": 126.904,  "group": 17, "period": 5, "block": "p",  "en": 2.66, "config": "[Kr] 4d10 5s2 5p5", "category": "halogen"},
    "Xe": {"Z": 54, "name": "Xenon",        "mass": 131.293,  "group": 18, "period": 5, "block": "p",  "en": 2.60, "config": "[Kr] 4d10 5s2 5p6", "category": "noble gas"},
    "Cs": {"Z": 55, "name": "Cesium",       "mass": 132.905,  "group": 1,  "period": 6, "block": "s",  "en": 0.79, "config": "[Xe] 6s1",      "category": "alkali metal"},
    "Ba": {"Z": 56, "name": "Barium",       "mass": 137.327,  "group": 2,  "period": 6, "block": "s",  "en": 0.89, "config": "[Xe] 6s2",      "category": "alkaline earth"},
    "W":  {"Z": 74, "name": "Tungsten",     "mass": 183.84,   "group": 6,  "period": 6, "block": "d",  "en": 2.36, "config": "[Xe] 4f14 5d4 6s2", "category": "transition metal"},
    "Re": {"Z": 75, "name": "Rhenium",      "mass": 186.207,  "group": 7,  "period": 6, "block": "d",  "en": 1.90, "config": "[Xe] 4f14 5d5 6s2", "category": "transition metal"},
    "Os": {"Z": 76, "name": "Osmium",       "mass": 190.23,   "group": 8,  "period": 6, "block": "d",  "en": 2.20, "config": "[Xe] 4f14 5d6 6s2", "category": "transition metal"},
    "Ir": {"Z": 77, "name": "Iridium",      "mass": 192.217,  "group": 9,  "period": 6, "block": "d",  "en": 2.20, "config": "[Xe] 4f14 5d7 6s2", "category": "transition metal"},
    "Pt": {"Z": 78, "name": "Platinum",     "mass": 195.084,  "group": 10, "period": 6, "block": "d",  "en": 2.28, "config": "[Xe] 4f14 5d9 6s1", "category": "transition metal"},
    "Au": {"Z": 79, "name": "Gold",         "mass": 196.967,  "group": 11, "period": 6, "block": "d",  "en": 2.54, "config": "[Xe] 4f14 5d10 6s1", "category": "transition metal"},
    "Hg": {"Z": 80, "name": "Mercury",      "mass": 200.59,   "group": 12, "period": 6, "block": "d",  "en": 2.00, "config": "[Xe] 4f14 5d10 6s2", "category": "transition metal"},
    "Tl": {"Z": 81, "name": "Thallium",     "mass": 204.38,   "group": 13, "period": 6, "block": "p",  "en": 1.62, "config": "[Xe] 4f14 5d10 6s2 6p1", "category": "post-transition metal"},
    "Pb": {"Z": 82, "name": "Lead",         "mass": 207.2,    "group": 14, "period": 6, "block": "p",  "en": 2.33, "config": "[Xe] 4f14 5d10 6s2 6p2", "category": "post-transition metal"},
    "Bi": {"Z": 83, "name": "Bismuth",      "mass": 208.980,  "group": 15, "period": 6, "block": "p",  "en": 2.02, "config": "[Xe] 4f14 5d10 6s2 6p3", "category": "post-transition metal"},
    "U":  {"Z": 92, "name": "Uranium",      "mass": 238.029,  "group": None, "period": 7, "block": "f", "en": 1.38, "config": "[Rn] 5f3 6d1 7s2", "category": "actinide"},
}

# Build reverse lookup by atomic number
_Z_TO_SYMBOL = {v["Z"]: k for k, v in ELEMENTS.items()}


def lookup_element(query: str) -> dict | None:
    """
    Look up element by symbol, name, or atomic number.
    Returns full element data dict or None.
    """
    q = query.strip()
    # By atomic number
    if q.isdigit():
        sym = _Z_TO_SYMBOL.get(int(q))
        if sym:
            return {"symbol": sym, **ELEMENTS[sym]}
        return None
    # By symbol (case-insensitive, capitalize first letter)
    sym = q.capitalize() if len(q) <= 2 else q
    if sym in ELEMENTS:
        return {"symbol": sym, **ELEMENTS[sym]}
    # By name
    q_lower = q.lower()
    for sym, data in ELEMENTS.items():
        if data["name"].lower() == q_lower:
            return {"symbol": sym, **data}
    # Fuzzy match
    matches = []
    for sym, data in ELEMENTS.items():
        if q_lower in data["name"].lower() or q_lower in sym.lower():
            matches.append({"symbol": sym, **data})
    return matches[0] if len(matches) == 1 else (matches if matches else None)


# =============================================================================
# CAS registry number validation
# =============================================================================


def validate_cas(cas_string: str) -> dict:
    """
    Validate a CAS registry number (format: NNNNNNN-NN-N).

    Returns dict with 'valid', 'cas', 'error' keys.
    """
    cas_string = cas_string.strip().replace(" ", "")

    # Accept with or without hyphens
    if "-" in cas_string:
        parts = cas_string.split("-")
        if len(parts) != 3:
            return {"valid": False, "cas": cas_string, "error": "CAS must have format NNNNN-NN-N"}
        try:
            part1, part2, check = parts
            if not (part1.isdigit() and part2.isdigit() and check.isdigit()):
                raise ValueError
            if len(part2) != 2 or len(check) != 1:
                raise ValueError
            digits = part1 + part2
            check_digit = int(check)
        except (ValueError, IndexError):
            return {"valid": False, "cas": cas_string, "error": "Invalid CAS format"}
    else:
        if not cas_string.isdigit() or len(cas_string) < 5:
            return {"valid": False, "cas": cas_string, "error": "Invalid CAS format"}
        digits = cas_string[:-1]
        check_digit = int(cas_string[-1])
        # Reconstruct hyphenated form
        cas_string = f"{digits[:-2]}-{digits[-2:]}-{check_digit}"

    # Validate check digit: sum of (position * digit) mod 10, counting from right
    total = sum((i + 1) * int(d) for i, d in enumerate(reversed(digits)))
    computed_check = total % 10

    if computed_check != check_digit:
        return {
            "valid": False,
            "cas": cas_string,
            "error": f"Check digit mismatch: expected {check_digit}, computed {computed_check}",
        }

    return {"valid": True, "cas": cas_string, "error": None}


# =============================================================================
# Isotope pattern calculation
# =============================================================================

# Natural isotope abundances (mass, abundance) for common elements
_ISOTOPES: dict[str, list[tuple[float, float]]] = {
    "H":  [(1.00783, 0.999885), (2.01410, 0.000115)],
    "C":  [(12.00000, 0.9893),  (13.00335, 0.0107)],
    "N":  [(14.00307, 0.99632), (15.00011, 0.00368)],
    "O":  [(15.99491, 0.99757), (16.99913, 0.00038), (17.99916, 0.00205)],
    "S":  [(31.97207, 0.9493),  (32.97146, 0.0076),  (33.96787, 0.0429), (35.96708, 0.0002)],
    "P":  [(30.97376, 1.0)],
    "F":  [(18.99840, 1.0)],
    "Cl": [(34.96885, 0.7576),  (36.96590, 0.2424)],
    "Br": [(78.91834, 0.5069),  (80.91629, 0.4931)],
    "I":  [(126.90447, 1.0)],
    "Si": [(27.97693, 0.92223), (28.97649, 0.04685), (29.97377, 0.03092)],
    "Se": [(73.92248, 0.0089),  (75.91921, 0.0937),  (76.91991, 0.0763),
           (77.91731, 0.2377),  (79.91652, 0.4961),  (81.91670, 0.0873)],
    "Na": [(22.98977, 1.0)],
    "K":  [(38.96371, 0.932581), (39.96400, 0.000117), (40.96183, 0.067302)],
    "Fe": [(53.93961, 0.05845), (55.93494, 0.91754), (56.93540, 0.02119), (57.93328, 0.00282)],
    "Cu": [(62.92960, 0.6915), (64.92779, 0.3085)],
    "Zn": [(63.92914, 0.4863), (65.92603, 0.2790), (66.92713, 0.0410),
           (67.92485, 0.1875), (69.92532, 0.0062)],
}

# Regex for molecular formula parsing: e.g., C14H19N3O4, Ca(OH)2, etc.
_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def _parse_formula_to_elements(formula: str) -> dict[str, int]:
    """Parse molecular formula string → {element: count}."""
    elements: dict[str, int] = {}
    # Simple parser (handles flat formulas like C9H8O4, not nested parentheses)
    # For parentheses, expand first
    expanded = _expand_parentheses(formula)
    for match in _FORMULA_RE.finditer(expanded):
        elem = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        if elem:
            elements[elem] = elements.get(elem, 0) + count
    return elements


def _expand_parentheses(formula: str) -> str:
    """Expand parenthesized groups: Ca(OH)2 → CaO2H2."""
    while "(" in formula:
        # Find innermost parentheses
        m = re.search(r"\(([^()]+)\)(\d*)", formula)
        if not m:
            break
        inner = m.group(1)
        mult = int(m.group(2)) if m.group(2) else 1
        # Multiply each element count
        expanded = ""
        for em in _FORMULA_RE.finditer(inner):
            elem = em.group(1)
            count = int(em.group(2)) if em.group(2) else 1
            if elem:
                expanded += f"{elem}{count * mult}"
        formula = formula[:m.start()] + expanded + formula[m.end():]
    return formula


def calculate_isotope_pattern(
    formula: str = "",
    smiles: str = "",
    charge: int = 1,
    top_n: int = 10,
    min_abundance: float = 0.001,
) -> dict:
    """
    Calculate isotope distribution from molecular formula or SMILES.

    Args:
        formula: Molecular formula (e.g., 'C9H8O4')
        smiles: SMILES string (alternative to formula, uses RDKit)
        charge: Charge state for m/z calculation (default 1)
        top_n: Max number of peaks to return
        min_abundance: Minimum relative abundance to include

    Returns dict with monoisotopic_mass, average_mass, pattern (list of {mz, abundance}).
    """
    if smiles and not formula:
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": f"Invalid SMILES: {smiles}"}
            formula = rdMolDescriptors.CalcMolFormula(mol)
        except ImportError:
            return {"error": "RDKit not available; provide formula instead"}

    if not formula:
        return {"error": "Provide either formula or smiles"}

    elements = _parse_formula_to_elements(formula)
    if not elements:
        return {"error": f"Could not parse formula: {formula}"}

    # Calculate pattern using polynomial multiplication
    pattern = [(0.0, 1.0)]  # (mass_offset, probability)

    for elem, count in elements.items():
        isotopes = _ISOTOPES.get(elem)
        if not isotopes:
            # Monoisotopic only for unknown elements
            mass = ELEMENTS.get(elem, {}).get("mass", 0)
            isotopes = [(mass, 1.0)]

        # For each atom of this element, convolve with isotope distribution
        for _ in range(count):
            new_pattern: dict[float, float] = {}
            for mass_a, prob_a in pattern:
                for mass_iso, prob_iso in isotopes:
                    combined_mass = round(mass_a + mass_iso, 5)
                    new_pattern[combined_mass] = new_pattern.get(combined_mass, 0.0) + prob_a * prob_iso
            # Prune low-abundance peaks for efficiency
            pattern = [(m, p) for m, p in new_pattern.items() if p > 1e-8]

    # Sort by mass, normalize
    pattern.sort(key=lambda x: x[0])
    max_prob = max(p for _, p in pattern) if pattern else 1.0

    # Aggregate peaks within 0.01 Da (numerical noise)
    aggregated: list[tuple[float, float]] = []
    for mass, prob in pattern:
        if aggregated and abs(mass - aggregated[-1][0]) < 0.01:
            # Merge
            old_m, old_p = aggregated[-1]
            total_p = old_p + prob
            avg_m = (old_m * old_p + mass * prob) / total_p
            aggregated[-1] = (avg_m, total_p)
        else:
            aggregated.append((mass, prob))

    max_prob = max(p for _, p in aggregated) if aggregated else 1.0

    result_peaks = []
    for mass, prob in aggregated:
        rel_abundance = prob / max_prob
        if rel_abundance >= min_abundance and len(result_peaks) < top_n:
            mz = mass / abs(charge) if charge != 0 else mass
            result_peaks.append({
                "mz": round(mz, 4),
                "abundance": round(prob * 100, 4),
                "relative": round(rel_abundance * 100, 2),
            })

    mono_mass = aggregated[0][0] if aggregated else 0.0
    avg_mass = sum(m * p for m, p in aggregated) / sum(p for _, p in aggregated) if aggregated else 0.0

    return {
        "formula": formula,
        "charge": charge,
        "monoisotopic_mass": round(mono_mass, 4),
        "monoisotopic_mz": round(mono_mass / abs(charge), 4) if charge != 0 else round(mono_mass, 4),
        "average_mass": round(avg_mass, 4),
        "pattern": result_peaks,
    }


# =============================================================================
# Scientific unit conversion
# =============================================================================

# Conversion factors to SI base units
_MASS_TO_KG = {
    "kg": 1, "g": 1e-3, "mg": 1e-6, "μg": 1e-9, "ug": 1e-9, "ng": 1e-12,
    "pg": 1e-15, "fg": 1e-18, "lb": 0.453592, "oz": 0.0283495,
    "Da": 1.66054e-27, "kDa": 1.66054e-24, "amu": 1.66054e-27,
}
_VOLUME_TO_L = {
    "L": 1, "l": 1, "mL": 1e-3, "ml": 1e-3, "μL": 1e-6, "uL": 1e-6,
    "nL": 1e-9, "pL": 1e-12, "dL": 0.1, "kL": 1000,
    "gal": 3.78541, "qt": 0.946353, "pt": 0.473176, "fl_oz": 0.0295735,
    "cm3": 1e-3, "m3": 1000, "mm3": 1e-6,
}
_LENGTH_TO_M = {
    "m": 1, "km": 1000, "cm": 0.01, "mm": 1e-3, "μm": 1e-6, "um": 1e-6,
    "nm": 1e-9, "pm": 1e-12, "Å": 1e-10, "A": 1e-10,
    "in": 0.0254, "ft": 0.3048, "mi": 1609.34,
}
_ENERGY_TO_J = {
    "J": 1, "kJ": 1000, "cal": 4.184, "kcal": 4184, "eV": 1.60218e-19,
    "keV": 1.60218e-16, "MeV": 1.60218e-13,
    "kJ/mol": 1000/6.02214e23, "kcal/mol": 4184/6.02214e23,
    "Eh": 4.35975e-18, "hartree": 4.35975e-18,
    "cm-1": 1.98645e-23, "wavenumber": 1.98645e-23,
}
_PRESSURE_TO_PA = {
    "Pa": 1, "kPa": 1000, "MPa": 1e6, "bar": 1e5, "mbar": 100,
    "atm": 101325, "torr": 133.322, "mmHg": 133.322, "psi": 6894.76,
}
_TIME_TO_S = {
    "s": 1, "ms": 1e-3, "μs": 1e-6, "us": 1e-6, "ns": 1e-9, "ps": 1e-12,
    "fs": 1e-15, "min": 60, "h": 3600, "hr": 3600, "d": 86400, "day": 86400,
}
_AMOUNT_TO_MOL = {
    "mol": 1, "mmol": 1e-3, "μmol": 1e-6, "umol": 1e-6,
    "nmol": 1e-9, "pmol": 1e-12, "fmol": 1e-15, "kmol": 1000,
}

_UNIT_CATEGORIES = {
    "mass": _MASS_TO_KG,
    "volume": _VOLUME_TO_L,
    "length": _LENGTH_TO_M,
    "energy": _ENERGY_TO_J,
    "pressure": _PRESSURE_TO_PA,
    "time": _TIME_TO_S,
    "amount": _AMOUNT_TO_MOL,
}


def _find_category(unit: str) -> tuple[str, dict] | None:
    """Find which category a unit belongs to."""
    for cat, table in _UNIT_CATEGORIES.items():
        if unit in table:
            return cat, table
    return None


def convert_units(value: float, from_unit: str, to_unit: str) -> dict:
    """
    Convert between scientific units.

    Supports: mass (g↔kg↔mg↔Da), volume (L↔mL↔μL), length (m↔nm↔Å),
    energy (J↔kJ↔kcal↔eV↔cm⁻¹↔hartree), pressure (Pa↔atm↔bar↔torr),
    time (s↔ms↔min↔h), amount (mol↔mmol↔μmol), and temperature (°C↔°F↔K).
    """
    # Temperature special case
    temp_units = {"C", "°C", "F", "°F", "K"}
    fu = from_unit.strip()
    tu = to_unit.strip()

    if fu in temp_units or tu in temp_units:
        fu_norm = fu.replace("°", "")
        tu_norm = tu.replace("°", "")
        # Convert to Kelvin first
        if fu_norm == "C":
            k = value + 273.15
        elif fu_norm == "F":
            k = (value - 32) * 5/9 + 273.15
        elif fu_norm == "K":
            k = value
        else:
            return {"error": f"Unknown temperature unit: {fu}"}
        # Convert from Kelvin
        if tu_norm == "C":
            result = k - 273.15
        elif tu_norm == "F":
            result = (k - 273.15) * 9/5 + 32
        elif tu_norm == "K":
            result = k
        else:
            return {"error": f"Unknown temperature unit: {tu}"}
        return {
            "value": value, "from_unit": from_unit,
            "result": round(result, 6), "to_unit": to_unit,
            "category": "temperature",
        }

    cat_from = _find_category(fu)
    cat_to = _find_category(tu)

    if not cat_from:
        return {"error": f"Unknown unit: {fu}. Supported: {', '.join(u for table in _UNIT_CATEGORIES.values() for u in table)}"}
    if not cat_to:
        return {"error": f"Unknown unit: {tu}. Supported: {', '.join(u for table in _UNIT_CATEGORIES.values() for u in table)}"}

    if cat_from[0] != cat_to[0]:
        return {"error": f"Cannot convert between {cat_from[0]} ({fu}) and {cat_to[0]} ({tu})"}

    # Convert: value * (from_factor / to_factor)
    si_value = value * cat_from[1][fu]
    result = si_value / cat_to[1][tu]

    return {
        "value": value, "from_unit": from_unit,
        "result": result, "to_unit": to_unit,
        "category": cat_from[0],
    }


# =============================================================================
# Buffer pH / Henderson-Hasselbalch
# =============================================================================

# Common buffer pKa values at 25°C
BUFFER_PKA: dict[str, dict] = {
    "phosphate_1":   {"name": "Phosphoric acid (pKa₁)", "pKa": 2.15, "species": "H₃PO₄ / H₂PO₄⁻"},
    "citrate_1":     {"name": "Citric acid (pKa₁)",     "pKa": 3.13, "species": "H₃Cit / H₂Cit⁻"},
    "formate":       {"name": "Formic acid",             "pKa": 3.75, "species": "HCOOH / HCOO⁻"},
    "citrate_2":     {"name": "Citric acid (pKa₂)",     "pKa": 4.76, "species": "H₂Cit⁻ / HCit²⁻"},
    "acetate":       {"name": "Acetic acid",             "pKa": 4.76, "species": "CH₃COOH / CH₃COO⁻"},
    "citrate_3":     {"name": "Citric acid (pKa₃)",     "pKa": 6.40, "species": "HCit²⁻ / Cit³⁻"},
    "MES":           {"name": "MES",                     "pKa": 6.15, "species": "MES-H / MES"},
    "PIPES":         {"name": "PIPES",                   "pKa": 6.76, "species": "PIPES-H / PIPES"},
    "phosphate_2":   {"name": "Phosphoric acid (pKa₂)", "pKa": 7.20, "species": "H₂PO₄⁻ / HPO₄²⁻"},
    "MOPS":          {"name": "MOPS",                    "pKa": 7.20, "species": "MOPS-H / MOPS"},
    "HEPES":         {"name": "HEPES",                   "pKa": 7.55, "species": "HEPES-H / HEPES"},
    "imidazole":     {"name": "Imidazole",               "pKa": 6.99, "species": "ImH⁺ / Im"},
    "Tris":          {"name": "Tris",                    "pKa": 8.07, "species": "TrisH⁺ / Tris"},
    "TAPS":          {"name": "TAPS",                    "pKa": 8.44, "species": "TAPS-H / TAPS"},
    "borate":        {"name": "Boric acid",              "pKa": 9.24, "species": "B(OH)₃ / B(OH)₄⁻"},
    "CHES":          {"name": "CHES",                    "pKa": 9.50, "species": "CHES-H / CHES"},
    "glycine_amino": {"name": "Glycine (amino)",         "pKa": 9.60, "species": "Gly⁺ / Gly⁻"},
    "CAPS":          {"name": "CAPS",                    "pKa": 10.40, "species": "CAPS-H / CAPS"},
    "carbonate_1":   {"name": "Carbonic acid (pKa₁)",   "pKa": 6.35, "species": "H₂CO₃ / HCO₃⁻"},
    "carbonate_2":   {"name": "Carbonic acid (pKa₂)",   "pKa": 10.33, "species": "HCO₃⁻ / CO₃²⁻"},
    "phosphate_3":   {"name": "Phosphoric acid (pKa₃)", "pKa": 12.32, "species": "HPO₄²⁻ / PO₄³⁻"},
}


def calculate_buffer_ph(
    pKa: float | None = None,
    buffer_name: str | None = None,
    acid_conc: float | None = None,
    base_conc: float | None = None,
    ratio_base_acid: float | None = None,
    target_ph: float | None = None,
) -> dict:
    """
    Henderson-Hasselbalch calculations.

    Mode 1: Given pKa + concentrations/ratio → calculate pH
    Mode 2: Given pKa + target pH → calculate required ratio

    Args:
        pKa: pKa of the buffer (or provide buffer_name to look up)
        buffer_name: Look up pKa from known buffers (e.g., 'Tris', 'HEPES', 'phosphate_2')
        acid_conc: Concentration of acid form (any consistent unit)
        base_conc: Concentration of base form
        ratio_base_acid: [A⁻]/[HA] ratio (alternative to concentrations)
        target_ph: Target pH to calculate required ratio
    """
    # Resolve pKa
    if buffer_name and pKa is None:
        buf = BUFFER_PKA.get(buffer_name)
        if not buf:
            # Fuzzy search
            q = buffer_name.lower()
            matches = [(k, v) for k, v in BUFFER_PKA.items() if q in k.lower() or q in v["name"].lower()]
            if matches:
                buf = matches[0][1]
            else:
                avail = ", ".join(BUFFER_PKA.keys())
                return {"error": f"Buffer '{buffer_name}' not found. Available: {avail}"}
        pKa = buf["pKa"]
        species = buf.get("species", "")
    else:
        species = ""

    if pKa is None:
        return {"error": "Provide pKa or buffer_name"}

    result: dict[str, Any] = {
        "pKa": pKa,
        "buffer": buffer_name,
        "species": species,
        "buffer_range": f"{pKa - 1:.1f} – {pKa + 1:.1f}",
    }

    # Mode 2: target pH → ratio
    if target_ph is not None:
        ratio = 10 ** (target_ph - pKa)
        result["target_pH"] = target_ph
        result["required_ratio_base_acid"] = round(ratio, 4)
        result["note"] = f"Mix [A⁻]/[HA] = {ratio:.4f} to achieve pH {target_ph}"
        if abs(target_ph - pKa) > 1:
            result["warning"] = f"Target pH {target_ph} is outside buffer range ({pKa-1:.1f}–{pKa+1:.1f}). Buffer capacity will be poor."
        return result

    # Mode 1: concentrations/ratio → pH
    if ratio_base_acid is None:
        if acid_conc is not None and base_conc is not None and acid_conc > 0:
            ratio_base_acid = base_conc / acid_conc
        else:
            return {"error": "Provide acid_conc + base_conc, or ratio_base_acid, or target_ph"}

    if ratio_base_acid <= 0:
        return {"error": "Ratio [base]/[acid] must be > 0"}

    ph = pKa + math.log10(ratio_base_acid)
    result["ratio_base_acid"] = round(ratio_base_acid, 4)
    result["calculated_pH"] = round(ph, 2)
    if acid_conc is not None:
        result["acid_concentration"] = acid_conc
    if base_conc is not None:
        result["base_concentration"] = base_conc
    if abs(ph - pKa) > 1:
        result["warning"] = f"pH {ph:.2f} is outside optimal buffer range ({pKa-1:.1f}–{pKa+1:.1f})"

    return result
