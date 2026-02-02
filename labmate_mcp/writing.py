"""
labmate-mcp writing & publication support module.

Tools for the final phase of the research workflow:
  - Citation formatting via Crossref content negotiation
  - Experimental section templates for common reaction types
  - Journal submission formatting guides
  - Supporting Information (SI) checklists
  - Standard chemistry abbreviations & symbols
  - Thesis writing section guides
  - Molecular formula formatting (LaTeX, HTML, Unicode)
  - IUPAC nomenclature lookup (PubChem)

No API keys required for most functions. Crossref and PubChem calls are free.
"""

from __future__ import annotations

import logging
import re
from typing import Any

import httpx

logger = logging.getLogger(__name__)


# =============================================================================
# Citation formatting via Crossref content negotiation
# =============================================================================

# CSL style → Crossref content-negotiation style parameter
CSL_STYLES: dict[str, str] = {
    "acs":          "american-chemical-society",
    "acs-nano":     "acs-nano",
    "jacs":         "american-chemical-society",
    "apa":          "apa",
    "apa7":         "apa-7th-edition",
    "rsc":          "royal-society-of-chemistry",
    "nature":       "nature",
    "science":      "science",
    "angew":        "angewandte-chemie",
    "angewandte":   "angewandte-chemie",
    "ieee":         "ieee",
    "vancouver":    "vancouver",
    "chicago":      "chicago-author-date",
    "harvard":      "harvard-cite-them-right",
    "cell":         "cell",
    "pnas":         "proceedings-of-the-national-academy-of-sciences",
    "elsevier":     "elsevier-harvard",
    "springer":     "springer-basic-author-date",
    "wiley":        "american-chemical-society",
    "mla":          "modern-language-association",
}


async def format_citation(
    doi: str,
    style: str = "acs",
    locale: str = "en-US",
) -> dict:
    """
    Format a DOI as a citation in the requested style using Crossref content negotiation.

    Args:
        doi: DOI string (with or without https://doi.org/ prefix)
        style: Citation style (acs, apa, rsc, nature, angew, ieee, vancouver, etc.)
        locale: Locale for formatting (default: en-US)

    Returns:
        dict with 'citation', 'style', 'doi' keys.
    """
    doi = doi.strip()
    if doi.startswith("http"):
        doi = doi.split("doi.org/")[-1]

    csl_style = CSL_STYLES.get(style.lower(), style.lower())

    url = f"https://api.crossref.org/works/{doi}/transform"
    headers = {
        "Accept": f"text/x-bibliography; style={csl_style}; locale={locale}",
        "User-Agent": "labmate-mcp/7.0.0 (mailto:labmate@scholarly.dev)",
    }

    try:
        async with httpx.AsyncClient(timeout=15, follow_redirects=True) as client:
            resp = await client.get(url, headers=headers)
            if resp.status_code == 200:
                citation = resp.text.strip()
                return {
                    "doi": doi,
                    "style": style,
                    "csl_style": csl_style,
                    "citation": citation,
                }
            elif resp.status_code == 404:
                return {"error": f"DOI not found: {doi}"}
            else:
                return {"error": f"Crossref returned {resp.status_code}", "doi": doi}
    except Exception as e:
        return {"error": f"Network error: {e}", "doi": doi}


async def build_bibliography(
    dois: list[str],
    style: str = "acs",
    numbered: bool = True,
    locale: str = "en-US",
) -> dict:
    """
    Build a formatted bibliography from a list of DOIs.

    Returns dict with 'bibliography' (formatted string) and 'entries' (individual citations).
    """
    entries = []
    errors = []

    for i, doi in enumerate(dois[:100]):  # cap at 100
        result = await format_citation(doi, style=style, locale=locale)
        if "error" in result:
            errors.append({"doi": doi, "error": result["error"]})
        else:
            entries.append({
                "index": i + 1,
                "doi": doi,
                "citation": result["citation"],
            })

    # Build formatted bibliography
    lines = []
    for entry in entries:
        if numbered:
            lines.append(f"({entry['index']}) {entry['citation']}")
        else:
            lines.append(entry["citation"])

    return {
        "style": style,
        "num_entries": len(entries),
        "num_errors": len(errors),
        "bibliography": "\n\n".join(lines),
        "entries": entries,
        "errors": errors if errors else None,
    }


# =============================================================================
# IUPAC name ↔ SMILES via PubChem PUG-REST
# =============================================================================

PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


async def iupac_from_smiles(smiles: str) -> dict:
    """Look up IUPAC name from a SMILES string via PubChem."""
    url = f"{PUBCHEM_REST}/compound/smiles/property/IUPACName,MolecularFormula,MolecularWeight,IsomericSMILES/JSON"
    params = {"smiles": smiles}
    try:
        async with httpx.AsyncClient(timeout=15) as client:
            resp = await client.get(url, params=params)
            if resp.status_code == 200:
                data = resp.json()
                props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
                return {
                    "input_smiles": smiles,
                    "iupac_name": props.get("IUPACName"),
                    "formula": props.get("MolecularFormula"),
                    "molecular_weight": props.get("MolecularWeight"),
                    "canonical_smiles": props.get("IsomericSMILES"),
                }
            else:
                return {"error": f"PubChem returned {resp.status_code}. Compound may not be in database.", "smiles": smiles}
    except Exception as e:
        return {"error": str(e), "smiles": smiles}


async def smiles_from_name(name: str) -> dict:
    """Look up SMILES and properties from a chemical name via PubChem."""
    import urllib.parse
    encoded = urllib.parse.quote(name, safe="")
    url = f"{PUBCHEM_REST}/compound/name/{encoded}/property/IsomericSMILES,CanonicalSMILES,IUPACName,MolecularFormula,MolecularWeight,InChI,InChIKey/JSON"
    try:
        async with httpx.AsyncClient(timeout=15) as client:
            resp = await client.get(url)
            if resp.status_code == 200:
                data = resp.json()
                props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
                return {
                    "input_name": name,
                    "isomeric_smiles": props.get("IsomericSMILES"),
                    "canonical_smiles": props.get("CanonicalSMILES"),
                    "iupac_name": props.get("IUPACName"),
                    "formula": props.get("MolecularFormula"),
                    "molecular_weight": props.get("MolecularWeight"),
                    "inchi": props.get("InChI"),
                    "inchi_key": props.get("InChIKey"),
                }
            elif resp.status_code == 404:
                return {"error": f"Compound not found: '{name}'"}
            else:
                return {"error": f"PubChem returned {resp.status_code}", "name": name}
    except Exception as e:
        return {"error": str(e), "name": name}


# =============================================================================
# Molecular formula formatting
# =============================================================================


def format_molecular_formula(
    formula: str,
    output_format: str = "unicode",
) -> dict:
    """
    Format a molecular formula with proper subscripts.

    Args:
        formula: e.g., "C9H8O4", "Ca(OH)2", "Fe2O3"
        output_format: "unicode", "latex", "html", or "plain"

    Returns dict with formatted string.
    """
    f = formula.strip()

    if output_format == "latex":
        # C9H8O4 → C$_{9}$H$_{8}$O$_{4}$
        out = re.sub(r"(\d+)", r"$_{\1}$", f)
        # Also handle in \ce{} notation
        ce = re.sub(r"(\d+)", r"_{\1}", f)
        return {"formula": f, "latex_inline": out, "latex_ce": f"\\ce{{{ce}}}", "format": "latex"}

    elif output_format == "html":
        out = re.sub(r"(\d+)", r"<sub>\1</sub>", f)
        return {"formula": f, "html": out, "format": "html"}

    elif output_format == "unicode":
        subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
        out = ""
        for ch in f:
            if ch.isdigit():
                out += ch.translate(subscript_map)
            else:
                out += ch
        return {"formula": f, "unicode": out, "format": "unicode"}

    else:  # plain
        return {"formula": f, "plain": f, "format": "plain"}


# =============================================================================
# Experimental section templates
# =============================================================================

EXPERIMENTAL_TEMPLATES: dict[str, dict] = {}


def _et(name, *, aliases=None, category="", template="", notes="", safety=""):
    key = name.lower()
    EXPERIMENTAL_TEMPLATES[key] = {
        "name": name,
        "aliases": aliases or [],
        "category": category,
        "template": template,
        "notes": notes,
        "safety": safety,
    }
    for a in (aliases or []):
        EXPERIMENTAL_TEMPLATES[a.lower()] = EXPERIMENTAL_TEMPLATES[key]


# --- Cross-Coupling Reactions ---

_et("Suzuki Coupling",
    aliases=["suzuki-miyaura", "suzuki"],
    category="Cross-Coupling",
    template="""**{product_name}**

A flame-dried round-bottom flask equipped with a magnetic stir bar was charged with {aryl_halide} ({aryl_halide_mass} mg, {aryl_halide_mmol} mmol, {aryl_halide_equiv} equiv), {boronic_acid} ({boronic_acid_mass} mg, {boronic_acid_mmol} mmol, {boronic_acid_equiv} equiv), Pd(PPh₃)₄ ({cat_mass} mg, {cat_mmol} mmol, {cat_mol_percent} mol%), and {base} ({base_mass} mg, {base_mmol} mmol, {base_equiv} equiv). The flask was evacuated and backfilled with N₂ (3×). {solvent} ({solvent_volume} mL) was added, and the reaction mixture was stirred at {temperature} °C for {time} h. The reaction was monitored by TLC ({tlc_system}). Upon completion, the mixture was cooled to room temperature, diluted with EtOAc ({workup_volume} mL), and washed with H₂O ({workup_volume} mL) and brine ({workup_volume} mL). The organic layer was dried over Na₂SO₄, filtered, and concentrated under reduced pressure. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound as a {product_appearance} ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.
{mp_or_rotation}""",
    notes="Replace placeholders with actual values. For Pd(dppf)Cl₂ catalyst, use degassed dioxane/H₂O (4:1).",
    safety="Pd catalysts: avoid inhalation. Boronic acids: irritants. Use fume hood.")

_et("Sonogashira Coupling",
    aliases=["sonogashira"],
    category="Cross-Coupling",
    template="""**{product_name}**

A Schlenk flask was charged with {aryl_halide} ({aryl_halide_mass} mg, {aryl_halide_mmol} mmol, {aryl_halide_equiv} equiv), PdCl₂(PPh₃)₂ ({pd_mass} mg, {pd_mmol} mmol, {pd_mol_percent} mol%), and CuI ({cui_mass} mg, {cui_mmol} mmol, {cui_mol_percent} mol%). The flask was evacuated and backfilled with N₂ (3×). Degassed {solvent} ({solvent_volume} mL) and {amine_base} ({amine_volume} mL) were added via syringe. {alkyne} ({alkyne_mass} mg, {alkyne_mmol} mmol, {alkyne_equiv} equiv) was then added dropwise. The mixture was stirred at {temperature} °C for {time} h. The reaction was filtered through Celite, washed with EtOAc, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Use degassed solvents. CuI amount typically 5-10 mol%. Et₃N or iPr₂NH as base.",
    safety="CuI: toxic. Pd compounds: avoid inhalation.")

_et("Buchwald-Hartwig Amination",
    aliases=["buchwald-hartwig", "c-n coupling"],
    category="Cross-Coupling",
    template="""**{product_name}**

An oven-dried Schlenk tube was charged with Pd₂(dba)₃ ({pd_mass} mg, {pd_mmol} mmol, {pd_mol_percent} mol%), {ligand} ({lig_mass} mg, {lig_mmol} mmol, {lig_equiv} equiv), and {base} ({base_mass} mg, {base_mmol} mmol, {base_equiv} equiv). The tube was evacuated and backfilled with N₂ (3×). {aryl_halide} ({aryl_halide_mass} mg, {aryl_halide_mmol} mmol, {aryl_halide_equiv} equiv), {amine} ({amine_mass} mg, {amine_mmol} mmol, {amine_equiv} equiv), and {solvent} ({solvent_volume} mL) were added. The reaction was heated to {temperature} °C and stirred for {time} h. After cooling, the mixture was diluted with EtOAc, filtered through Celite, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Common ligands: BINAP, XPhos, SPhos, RuPhos, DavePhos, BrettPhos. Base: NaOtBu or Cs₂CO₃.",
    safety="Pd catalysts: avoid inhalation. Strong bases: moisture-sensitive. Use glovebox or Schlenk line.")

# --- Reductions ---

_et("Sodium Borohydride Reduction",
    aliases=["nabh4 reduction", "nabh4", "borohydride reduction"],
    category="Reduction",
    template="""**{product_name}**

To a solution of {substrate} ({substrate_mass} mg, {substrate_mmol} mmol, {substrate_equiv} equiv) in {solvent} ({solvent_volume} mL) at {temperature} °C was added NaBH₄ ({nabh4_mass} mg, {nabh4_mmol} mmol, {nabh4_equiv} equiv) portionwise over {addition_time} min. The reaction was stirred at {temperature} °C for {time} h. The reaction was quenched by careful addition of saturated NH₄Cl solution ({quench_volume} mL) at 0 °C and extracted with {extraction_solvent} (3 × {extraction_volume} mL). The combined organic layers were washed with brine, dried over Na₂SO₄, filtered, and concentrated under reduced pressure. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound as a {product_appearance} ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Quench carefully — H₂ evolution. For ketone selectivity over ester, use MeOH at −78 °C. CeCl₃ for Luche conditions.",
    safety="NaBH₄: flammable solid, water-reactive. Quench generates H₂ — open vessel, no sparks.")

_et("LiAlH4 Reduction",
    aliases=["lah reduction", "lialh4", "lithium aluminium hydride reduction"],
    category="Reduction",
    template="""**{product_name}**

To a suspension of LiAlH₄ ({lah_mass} mg, {lah_mmol} mmol, {lah_equiv} equiv) in dry {solvent} ({solvent_volume} mL) at 0 °C under N₂ was added a solution of {substrate} ({substrate_mass} mg, {substrate_mmol} mmol, {substrate_equiv} equiv) in dry {solvent} ({substrate_solvent_volume} mL) dropwise via addition funnel over {addition_time} min. The reaction was allowed to warm to {temperature} °C and stirred for {time} h. The reaction was cooled to 0 °C and carefully quenched by sequential addition of H₂O ({fieser_water} mL), 15% NaOH ({fieser_naoh} mL), and H₂O ({fieser_water2} mL) [Fieser workup]. The resulting white precipitate was filtered through Celite and washed with {solvent}. The filtrate was dried over Na₂SO₄, filtered, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Fieser workup: for each gram of LiAlH₄, add x mL H₂O, x mL 15% NaOH, 3x mL H₂O. Alternative: Rochelle's salt (Na/K tartrate) workup.",
    safety="LiAlH₄: pyrophoric, water-reactive. Anhydrous conditions mandatory. Quench under N₂ at 0 °C. Fire extinguisher on hand.")

_et("Catalytic Hydrogenation",
    aliases=["hydrogenation", "pd/c hydrogenation", "h2 reduction"],
    category="Reduction",
    template="""**{product_name}**

A round-bottom flask was charged with {substrate} ({substrate_mass} mg, {substrate_mmol} mmol), Pd/C (10 wt%, {cat_mass} mg, {cat_loading} mol% Pd), and {solvent} ({solvent_volume} mL). The flask was evacuated and backfilled with H₂ (balloon, 3×). The reaction was stirred at room temperature under H₂ atmosphere (1 atm, balloon) for {time} h. The reaction was filtered through a pad of Celite, washed with {solvent} ({wash_volume} mL), and concentrated under reduced pressure to afford the title compound as a {product_appearance} ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Pd/C: 10% w/w typical, can use 5-20%. Solvents: MeOH, EtOAc, EtOH. For debenzylation, use same conditions. Pd(OH)₂/C (Pearlman's catalyst) for stubborn substrates.",
    safety="Pd/C: pyrophoric when dry. Always add solvent first. Never expose dry Pd/C to H₂. Filter away from open flames.")

# --- Oxidations ---

_et("Dess-Martin Oxidation",
    aliases=["dmp oxidation", "dess-martin"],
    category="Oxidation",
    template="""**{product_name}**

To a solution of {substrate} ({substrate_mass} mg, {substrate_mmol} mmol, 1.0 equiv) in CH₂Cl₂ ({solvent_volume} mL) at room temperature was added Dess-Martin periodinane ({dmp_mass} mg, {dmp_mmol} mmol, {dmp_equiv} equiv). The reaction was stirred at room temperature for {time} h, then quenched by addition of a 1:1 mixture of saturated NaHCO₃ and saturated Na₂S₂O₃ ({quench_volume} mL). The layers were separated, and the aqueous layer was extracted with CH₂Cl₂ (3 × {extraction_volume} mL). The combined organic layers were washed with brine, dried over Na₂SO₄, filtered, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="1.1-1.5 equiv DMP. Na₂S₂O₃ reduces excess periodinane and iodinane byproducts. No over-oxidation to carboxylic acid.",
    safety="DMP: oxidizer, shock-sensitive when dry. Store cold. Waste: treat with Na₂S₂O₃ before disposal.")

_et("Swern Oxidation",
    aliases=["swern"],
    category="Oxidation",
    template="""**{product_name}**

To a solution of oxalyl chloride ({oxchloride_volume} μL, {oxchloride_mmol} mmol, {oxchloride_equiv} equiv) in CH₂Cl₂ ({solvent_volume} mL) at −78 °C was added DMSO ({dmso_volume} μL, {dmso_mmol} mmol, {dmso_equiv} equiv) dropwise. The mixture was stirred for 15 min, then a solution of {substrate} ({substrate_mass} mg, {substrate_mmol} mmol, 1.0 equiv) in CH₂Cl₂ ({substrate_solvent_volume} mL) was added dropwise. The reaction was stirred at −78 °C for {time} min, then Et₃N ({et3n_volume} mL, {et3n_mmol} mmol, {et3n_equiv} equiv) was added. The reaction was allowed to warm to room temperature over 1 h, then diluted with CH₂Cl₂ ({dilute_volume} mL) and washed with H₂O, 1 M HCl, saturated NaHCO₃, and brine. The organic layer was dried over Na₂SO₄, filtered, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Order of addition critical: (COCl)₂ first, then DMSO, then substrate, then Et₃N. Keep at −78 °C until Et₃N added.",
    safety="Oxalyl chloride: highly toxic, corrosive, lachrymator. DMSO/oxalyl chloride: exothermic at >−60 °C. Use dry ice/acetone bath. Fume hood mandatory.")

# --- Amide Coupling ---

_et("Amide Coupling",
    aliases=["peptide coupling", "hatu coupling", "edc coupling", "amidation"],
    category="Amide Bond Formation",
    template="""**{product_name}**

To a solution of {acid} ({acid_mass} mg, {acid_mmol} mmol, {acid_equiv} equiv) in {solvent} ({solvent_volume} mL) at {temperature} °C was added {coupling_reagent} ({coupling_mass} mg, {coupling_mmol} mmol, {coupling_equiv} equiv) and {base} ({base_volume} μL, {base_mmol} mmol, {base_equiv} equiv). The mixture was stirred for {preactivation_time} min, then {amine} ({amine_mass} mg, {amine_mmol} mmol, {amine_equiv} equiv) was added. The reaction was stirred at room temperature for {time} h. The reaction was diluted with EtOAc ({dilute_volume} mL) and washed sequentially with 1 M HCl, saturated NaHCO₃, and brine. The organic layer was dried over Na₂SO₄, filtered, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Coupling reagents: HATU > HBTU > TBTU > EDC·HCl/HOBt. Bases: DIPEA (iPr₂NEt), NMM, Et₃N. Solvents: DMF, CH₂Cl₂, or DMF/CH₂Cl₂ mixtures. Pre-activation with coupling reagent + base for 5-15 min improves yields.",
    safety="HATU/HBTU: skin sensitizers, irritants. DMF: reproductive toxicant — avoid skin contact.")

# --- SPPS ---

_et("Solid-Phase Peptide Synthesis",
    aliases=["spps", "fmoc spps", "peptide synthesis"],
    category="Peptide Chemistry",
    template="""**Solid-Phase Peptide Synthesis of {peptide_name}**

Peptide synthesis was performed manually on a {resin_type} resin (loading: {resin_loading} mmol/g, {resin_mass} mg, {scale} mmol scale) using standard Fmoc/tBu chemistry.

*Fmoc Deprotection*: The resin was treated with 20% piperidine in DMF (2 × {deprot_volume} mL, 5 + 15 min) and washed with DMF (5 × {wash_volume} mL).

*Amino Acid Coupling*: Fmoc-{aa_name}-OH ({aa_equiv} equiv), {coupling_reagent} ({coupling_equiv} equiv), and {base} ({base_equiv} equiv) were dissolved in DMF ({coupling_volume} mL) and added to the resin. The mixture was agitated for {coupling_time} min at room temperature. Coupling completion was monitored by the Kaiser (ninhydrin) test. If positive, the coupling was repeated.

*Wash Cycles*: Between each step: DMF (5×), CH₂Cl₂ (3×), DMF (3×).

*Cleavage*: The peptide was cleaved from the resin using TFA/{scavenger_cocktail} ({cleavage_ratio}, {cleavage_volume} mL) for {cleavage_time} h at room temperature. The resin was filtered off, and the filtrate was concentrated under N₂ flow. The crude peptide was precipitated with cold diethyl ether ({ether_volume} mL), centrifuged (4000 rpm, 5 min), and the pellet was washed with cold ether (2×). The crude peptide was dissolved in {dissolve_solvent} and lyophilized to afford a {product_appearance} ({crude_mass} mg).

*Purification*: The crude peptide was purified by preparative RP-HPLC (Column: {hplc_column}; gradient: {hplc_gradient}; flow rate: {hplc_flow} mL/min; detection: λ = 220/280 nm). Pure fractions were pooled and lyophilized to afford the title peptide as a {pure_appearance} ({pure_mass} mg, {yield_percent}% overall).

Analytical HPLC: tR = {hplc_rt} min (purity: {hplc_purity}%, {hplc_conditions}).
MALDI-TOF MS: m/z calcd for {ms_formula} [M+H]⁺ {ms_calcd}, found {ms_found}.
{additional_characterization}""",
    notes="Scavenger cocktails: TFA/TIS/H₂O (95:2.5:2.5) standard. For Cys-containing: TFA/TIS/H₂O/EDT (94:1:2.5:2.5). For Met/Trp-containing: add EDT or DODT. Kaiser test: ninhydrin/pyridine/KCN; blue = free amine (incomplete coupling), yellow = complete.",
    safety="TFA: highly corrosive, strong acid. Piperidine: highly toxic, flammable. EDT: foul-smelling thiol. Use fume hood for all operations.")

# --- Click Chemistry ---

_et("CuAAC Click Reaction",
    aliases=["click chemistry", "copper click", "cuaac", "azide-alkyne cycloaddition"],
    category="Cycloaddition",
    template="""**{product_name}**

To a solution of {azide} ({azide_mass} mg, {azide_mmol} mmol, {azide_equiv} equiv) and {alkyne} ({alkyne_mass} mg, {alkyne_mmol} mmol, {alkyne_equiv} equiv) in {solvent} ({solvent_volume} mL) was added CuSO₄·5H₂O ({cuso4_mass} mg, {cuso4_mmol} mmol, {cu_mol_percent} mol%) and sodium ascorbate ({asc_mass} mg, {asc_mmol} mmol, {asc_equiv} equiv). The mixture was stirred at room temperature for {time} h. The reaction was diluted with EtOAc ({dilute_volume} mL), washed with saturated NH₄Cl ({wash_volume} mL), H₂O, and brine. The organic layer was dried over Na₂SO₄, filtered, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the 1,4-disubstituted 1,2,3-triazole product ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="CuSO₄/sodium ascorbate (in situ Cu(I)). Typically 5-10 mol% Cu. Solvents: t-BuOH/H₂O (1:1), DMF/H₂O, DMSO. TBTA ligand accelerates. Exclusively 1,4-regioisomer.",
    safety="Organic azides: potentially explosive — never isolate low-MW azides. CuSO₄: irritant. Ascorbate: safe.")

# --- Protection/Deprotection ---

_et("Boc Deprotection",
    aliases=["boc removal", "tfa deprotection", "boc cleavage"],
    category="Protection/Deprotection",
    template="""**{product_name}**

To a solution of {substrate} ({substrate_mass} mg, {substrate_mmol} mmol) in CH₂Cl₂ ({solvent_volume} mL) at 0 °C was added TFA ({tfa_volume} mL, {tfa_equiv} equiv). The reaction was stirred at room temperature for {time} h. The reaction was concentrated under reduced pressure. The residue was azeotroped with toluene (3×) to remove residual TFA, affording the TFA salt of the title compound as a {product_appearance} ({product_mass} mg, {yield_percent}%).

{characterization}""",
    notes="Typical: 25-50% TFA in CH₂Cl₂ (v/v). For acid-sensitive substrates: 4 M HCl in dioxane or TMSOTf/2,6-lutidine. Free-base by dissolving in CH₂Cl₂ and washing with saturated NaHCO₃.",
    safety="TFA: highly corrosive. Use fume hood. Neutralize waste with NaHCO₃ before disposal.")

_et("TBS Protection",
    aliases=["tbs silylation", "tbdms protection", "silyl protection"],
    category="Protection/Deprotection",
    template="""**{product_name}**

To a solution of {substrate} ({substrate_mass} mg, {substrate_mmol} mmol, 1.0 equiv) and imidazole ({imid_mass} mg, {imid_mmol} mmol, {imid_equiv} equiv) in {solvent} ({solvent_volume} mL) at {temperature} °C was added TBSCl ({tbs_mass} mg, {tbs_mmol} mmol, {tbs_equiv} equiv). The reaction was stirred at room temperature for {time} h. The reaction was quenched with H₂O ({quench_volume} mL) and extracted with {extraction_solvent} (3 × {extraction_volume} mL). The combined organic layers were washed with brine, dried over Na₂SO₄, filtered, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="TBSCl/imidazole in DMF for 1° alcohols. For 2°: TBSOTf/2,6-lutidine in CH₂Cl₂ at −78°C to 0°C. DMAP catalytic amount can accelerate.",
    safety="TBSCl: corrosive. Imidazole: irritant.")

_et("TBAF Desilylation",
    aliases=["tbaf deprotection", "tbs removal", "silyl removal"],
    category="Protection/Deprotection",
    template="""**{product_name}**

To a solution of {substrate} ({substrate_mass} mg, {substrate_mmol} mmol, 1.0 equiv) in THF ({solvent_volume} mL) at {temperature} °C was added TBAF (1.0 M in THF, {tbaf_volume} mL, {tbaf_mmol} mmol, {tbaf_equiv} equiv). The reaction was stirred at room temperature for {time} h, then quenched with saturated NH₄Cl ({quench_volume} mL). The mixture was extracted with EtOAc (3 × {extraction_volume} mL), washed with brine, dried over Na₂SO₄, filtered, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="1.0-2.0 equiv TBAF. Add AcOH to buffer if substrate is base-sensitive. Alternative: HF·pyridine (Olah's reagent) for hindered silyl ethers.",
    safety="TBAF: corrosive, hygroscopic. HF·pyridine: extremely dangerous — HF burns are life-threatening.")

# --- Recrystallization & Purification ---

_et("Recrystallization",
    aliases=["recrystallization", "crystallization"],
    category="Purification",
    template="""**Recrystallization of {compound_name}**

Crude {compound_name} ({crude_mass} mg) was dissolved in minimum hot {solvent} (~{dissolve_temp} °C, {solvent_volume} mL). The solution was allowed to cool slowly to room temperature, then placed at {cool_temp} °C for {cool_time} h. The resulting crystals were collected by vacuum filtration, washed with cold {wash_solvent} ({wash_volume} mL), and dried under vacuum to afford pure {compound_name} as {crystal_appearance} ({product_mass} mg, {recovery_percent}% recovery).

mp: {melting_point} °C.
{additional_characterization}""",
    notes="Solvent pairs: EtOAc/hexanes, CH₂Cl₂/hexanes, MeOH/H₂O, acetone/hexanes, EtOH/H₂O. Seed crystals improve nucleation. Slow cooling = larger, purer crystals.",
    safety="Use appropriate solvent safety precautions. Hot solvent: burn risk.")

# --- General Procedures ---

_et("Wittig Reaction",
    aliases=["wittig olefination", "wittig"],
    category="C-C Bond Formation",
    template="""**{product_name}**

To a suspension of {phosphonium_salt} ({salt_mass} mg, {salt_mmol} mmol, {salt_equiv} equiv) in THF ({solvent_volume} mL) at {base_temp} °C was added {base} ({base_volume} mL/mg, {base_mmol} mmol, {base_equiv} equiv) dropwise. The resulting orange/red ylide solution was stirred for {ylide_time} min. A solution of {aldehyde} ({aldehyde_mass} mg, {aldehyde_mmol} mmol, 1.0 equiv) in THF ({aldehyde_solvent_volume} mL) was added dropwise at {addition_temp} °C. The reaction was stirred at {reaction_temp} °C for {time} h. The reaction was quenched with saturated NH₄Cl, extracted with EtOAc (3×), washed with brine, dried over Na₂SO₄, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%, E:Z = {ez_ratio}).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Non-stabilized ylides: preferentially Z-alkene. Stabilized ylides: preferentially E-alkene. Bases: n-BuLi, NaHMDS, KHMDS. For E-selective olefination, use HWE (Horner-Wadsworth-Emmons) instead.",
    safety="n-BuLi: pyrophoric. Use syringe techniques under inert atmosphere.")

_et("Grignard Reaction",
    aliases=["grignard addition", "grignard"],
    category="C-C Bond Formation",
    template="""**{product_name}**

Mg turnings ({mg_mass} mg, {mg_mmol} mmol, {mg_equiv} equiv) were flame-dried under vacuum in a round-bottom flask. Dry THF ({thf_volume} mL) and a crystal of I₂ were added. {halide} ({halide_volume} μL, {halide_mmol} mmol, {halide_equiv} equiv) in THF ({halide_thf_volume} mL) was added dropwise at a rate to maintain gentle reflux. After addition, the mixture was refluxed for {grignard_time} h to ensure complete consumption of Mg. The Grignard reagent was cooled to {addition_temp} °C, and a solution of {electrophile} ({electrophile_mass} mg, {electrophile_mmol} mmol, 1.0 equiv) in THF ({electrophile_thf_volume} mL) was added dropwise. The reaction was stirred at {reaction_temp} °C for {time} h, then quenched with saturated NH₄Cl at 0 °C. The mixture was extracted with EtOAc (3×), washed with brine, dried over Na₂SO₄, and concentrated. Purification by column chromatography (SiO₂, {column_eluent}) afforded the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Anhydrous conditions essential. Activation: I₂, 1,2-dibromoethane, or DIBAL-H. Solvents: THF or Et₂O. Titrate Grignard with menthol/1,10-phenanthroline before use.",
    safety="Grignard formation is exothermic — control addition rate. Mg turnings: flammable. Anhydrous ether/THF: fire hazard.")

_et("Mitsunobu Reaction",
    aliases=["mitsunobu"],
    category="Substitution",
    template="""**{product_name}**

To a solution of {alcohol} ({alcohol_mass} mg, {alcohol_mmol} mmol, 1.0 equiv), {nucleophile} ({nuc_mass} mg, {nuc_mmol} mmol, {nuc_equiv} equiv), and PPh₃ ({pph3_mass} mg, {pph3_mmol} mmol, {pph3_equiv} equiv) in THF ({solvent_volume} mL) at 0 °C was added DIAD ({diad_volume} μL, {diad_mmol} mmol, {diad_equiv} equiv) dropwise. The reaction was allowed to warm to room temperature and stirred for {time} h. The reaction was concentrated and purified directly by column chromatography (SiO₂, {column_eluent}) to afford the title compound ({product_mass} mg, {yield_percent}%).

¹H NMR ({nmr_freq} MHz, {nmr_solvent}): δ {nmr_data_1h}
¹³C NMR ({nmr_freq_13c} MHz, {nmr_solvent}): δ {nmr_data_13c}
HRMS ({ms_method}): m/z calcd for {ms_formula} [M+{ms_ion}]⁺ {ms_calcd}, found {ms_found}.""",
    notes="Nucleophile pKₐ < 11 required. DIAD or DEAD as azodicarboxylate. Inversion of configuration at stereocenter. TPPO byproduct can be problematic — use polymer-bound PPh₃ or Fluorous-PPh₃ for easier removal.",
    safety="DIAD/DEAD: toxic, eye/skin irritant. PPh₃: sensitizer.")


def lookup_experimental_template(query: str) -> dict | list[dict] | None:
    """
    Look up an experimental section template by reaction name.

    Returns template dict or list of partial matches.
    """
    q = query.strip().lower()
    if q in EXPERIMENTAL_TEMPLATES:
        return EXPERIMENTAL_TEMPLATES[q]

    # Fuzzy search
    matches = []
    for key, tmpl in EXPERIMENTAL_TEMPLATES.items():
        if q in key or q in tmpl.get("category", "").lower():
            if tmpl not in matches:
                matches.append(tmpl)
    return matches if matches else None


# =============================================================================
# Journal formatting guides
# =============================================================================

JOURNAL_GUIDES: dict[str, dict] = {}


def _jg(name, *, aliases=None, publisher="", issn="", citation_style="",
        word_limit="", abstract_limit="", figure_limit="",
        reference_format="", si_allowed=True, graphical_abstract="",
        submission_url="", open_access="", apc="", scope="",
        file_formats="", special_notes=""):
    key = name.lower()
    JOURNAL_GUIDES[key] = {
        "name": name,
        "aliases": aliases or [],
        "publisher": publisher,
        "issn": issn,
        "citation_style": citation_style,
        "word_limit": word_limit,
        "abstract_limit": abstract_limit,
        "figure_limit": figure_limit,
        "reference_format": reference_format,
        "si_allowed": si_allowed,
        "graphical_abstract": graphical_abstract,
        "submission_url": submission_url,
        "open_access": open_access,
        "apc": apc,
        "scope": scope,
        "file_formats": file_formats,
        "special_notes": special_notes,
    }
    for a in (aliases or []):
        JOURNAL_GUIDES[a.lower()] = JOURNAL_GUIDES[key]


_jg("Journal of the American Chemical Society",
    aliases=["jacs", "j. am. chem. soc."],
    publisher="ACS", issn="0002-7863",
    citation_style="ACS (superscript numbered)",
    word_limit="Articles: no strict limit (typically 8000-12000 words). Communications: 4 journal pages (~3500 words + figures).",
    abstract_limit="200 words max",
    figure_limit="No strict limit for articles. Communications: typically 3-4 figures.",
    reference_format="ACS style: superscript numbers, numbered in order of appearance. Format: Author(s). Title. J. Abbrev. Year, Vol, Pages. DOI.",
    graphical_abstract="Required. 3.25 × 1.75 inches, max 600 dpi.",
    submission_url="https://pubs.acs.org/journal/jacsat",
    open_access="Hybrid. ACS AuthorChoice (CC-BY or CC-BY-NC-ND).",
    apc="~$5000 (ACS AuthorChoice)",
    scope="All areas of chemistry — emphasis on significance and broad interest.",
    file_formats="Manuscript: Word or LaTeX. Figures: TIFF, EPS, PDF (300+ dpi). SI: single PDF.",
    special_notes="TOC graphic required. ORCID mandatory for corresponding author. Data availability statement required.")

_jg("Angewandte Chemie International Edition",
    aliases=["angew", "angew. chem. int. ed.", "angewandte"],
    publisher="Wiley-VCH", issn="1433-7851",
    citation_style="Angew style: superscript numbered, numbered in order of citation.",
    word_limit="Communications: 4 pages (~2500 words + refs). Reviews: by invitation.",
    abstract_limit="200 words max",
    figure_limit="Communications: 3-5 figures/schemes/tables combined.",
    reference_format="Angew format: [N] a) Author, Journal Year, Vol, Pages; b) ... Multiple sub-references with a), b), c).",
    graphical_abstract="Table of Contents graphic required. Max 5 × 5 cm.",
    submission_url="https://onlinelibrary.wiley.com/journal/15213773",
    open_access="Hybrid. OnlineOpen (CC-BY, CC-BY-NC).",
    apc="~€4500",
    scope="All areas — must be of very high importance and broad interest.",
    file_formats="Word preferred. Figures: TIFF, EPS (300+ dpi).",
    special_notes="Very selective (~15% acceptance). Highlight (5 keywords) and catch-phrase for TOC required.")

_jg("Nature Chemistry",
    aliases=["nat. chem.", "nature chem"],
    publisher="Nature Publishing Group / Springer Nature", issn="1755-4330",
    citation_style="Nature style: superscript numbered, numbered in order of citation.",
    word_limit="Articles: ~3000 words (main text, excluding methods, refs, figure legends). Letters: ~1500 words.",
    abstract_limit="150 words max (no references in abstract)",
    figure_limit="Articles: up to 6 display items (figures/tables). Letters: up to 4.",
    reference_format="Nature style: Author(s). Title. Journal Abbrev. Vol, Pages (Year). Max 30 refs for Letters, 50 for Articles.",
    graphical_abstract="Not required (Nature uses own design for TOC).",
    submission_url="https://www.nature.com/nchem/",
    open_access="Hybrid. Gold OA available.",
    apc="~€9500 (Gold OA)",
    scope="All areas of chemistry — emphasis on novelty and broad significance.",
    file_formats="Word or LaTeX. Figures: separate files, 300+ dpi.",
    special_notes="Extremely selective. Methods section separate from main text. Extended Data figures allowed (8 max). Cover letter essential.")

_jg("Chemical Science",
    aliases=["chem. sci.", "chem sci"],
    publisher="Royal Society of Chemistry (RSC)", issn="2041-6520",
    citation_style="RSC style: superscript numbered.",
    word_limit="Edge articles: ~5 journal pages. Full articles: ~8 journal pages.",
    abstract_limit="200 words max",
    figure_limit="Reasonable (typically 4-8 figures).",
    reference_format="RSC format: N. Author, J. Abbrev., Year, Vol, Pages. ESI references numbered separately.",
    graphical_abstract="Table of Contents entry: max 8 × 4 cm.",
    submission_url="https://www.rsc.org/journals-books-databases/about-journals/chemical-science/",
    open_access="Full Gold OA (CC-BY). No APC — funded by RSC.",
    apc="Free (no APC)",
    scope="Chemistry of exceptional significance. All areas.",
    file_formats="Word or LaTeX. Figures: TIFF, EPS.",
    special_notes="Flagship RSC journal. Free open access (unique!). ESI (Electronic Supplementary Information) strongly encouraged.")

_jg("Organic Letters",
    aliases=["org. lett.", "orglett"],
    publisher="ACS", issn="1523-7060",
    citation_style="ACS style: superscript numbered.",
    word_limit="Communications only: max 4 printed pages (~2500 words).",
    abstract_limit="150 words max",
    figure_limit="Typically 2-4 schemes + 1-2 figures/tables.",
    reference_format="ACS style: superscript numbers.",
    graphical_abstract="TOC graphic required.",
    submission_url="https://pubs.acs.org/journal/orlef7",
    open_access="Hybrid (ACS AuthorChoice).",
    apc="~$5000",
    scope="Organic and bioorganic chemistry — synthesis, mechanisms, theory, new reagents.",
    file_formats="Word or LaTeX. SI: single PDF.",
    special_notes="Short communications only. Highly competitive. Experimental + full characterization in SI.")

_jg("The Journal of Organic Chemistry",
    aliases=["joc", "j. org. chem."],
    publisher="ACS", issn="0022-3263",
    citation_style="ACS style: superscript numbered.",
    word_limit="Full articles: no strict limit (typically 6000-12000). Notes: 4 pages.",
    abstract_limit="200 words max",
    figure_limit="Reasonable for full articles.",
    reference_format="ACS style: superscript numbers.",
    graphical_abstract="TOC graphic required.",
    submission_url="https://pubs.acs.org/journal/joceah",
    open_access="Hybrid.",
    apc="~$5000",
    scope="Organic chemistry: synthesis, mechanisms, theory, natural products, methodology.",
    file_formats="Word or LaTeX.",
    special_notes="Full experimental details in main manuscript (not SI). Complete NMR data for all new compounds.")

_jg("ACS Catalysis",
    aliases=["acs catal.", "acs catalysis"],
    publisher="ACS", issn="2155-5435",
    citation_style="ACS style.",
    word_limit="Letters: 3000 words. Articles: no strict limit.",
    abstract_limit="200 words",
    figure_limit="Reasonable.",
    reference_format="ACS style.",
    graphical_abstract="TOC graphic required.",
    submission_url="https://pubs.acs.org/journal/accacs",
    open_access="Hybrid.",
    apc="~$5000",
    scope="Heterogeneous, homogeneous, bio-, and electro-catalysis. Mechanisms and applications.",
    file_formats="Word or LaTeX.",
    special_notes="Scope includes computational catalysis. Turnover numbers/frequencies expected.")

_jg("Journal of Medicinal Chemistry",
    aliases=["jmc", "j. med. chem."],
    publisher="ACS", issn="0022-2623",
    citation_style="ACS style.",
    word_limit="Articles: 7500-10000 words. Letters: 3000 words. Perspectives: 10000 words.",
    abstract_limit="250 words",
    figure_limit="Reasonable.",
    reference_format="ACS style.",
    graphical_abstract="TOC graphic required.",
    submission_url="https://pubs.acs.org/journal/jmcmar",
    open_access="Hybrid.",
    apc="~$5000",
    scope="Drug design, SAR, ADME, computational med chem, chemical biology.",
    file_formats="Word or LaTeX.",
    special_notes="Requires biological data + SAR analysis. SMILES strings for all compounds. Molecular formula strings for HRMS.")

_jg("Chemical Communications",
    aliases=["chemcomm", "chem. commun.", "chem commun"],
    publisher="RSC", issn="1359-7345",
    citation_style="RSC style.",
    word_limit="Communications: 3500 words max (including refs and captions).",
    abstract_limit="No formal abstract — use first paragraph as summary.",
    figure_limit="4 figures/tables max.",
    reference_format="RSC format.",
    graphical_abstract="TOC entry + graphical abstract preferred.",
    submission_url="https://www.rsc.org/journals-books-databases/about-journals/chemcomm/",
    open_access="Hybrid. Gold OA available.",
    apc="~£2000",
    scope="All areas of chemistry. Short, urgent communications.",
    file_formats="Word. Figures: TIFF.",
    special_notes="Very short format. All experimental details must go in ESI. Rapid publication.")

_jg("Dalton Transactions",
    aliases=["dalton", "dalton trans."],
    publisher="RSC", issn="1477-9226",
    citation_style="RSC style.",
    word_limit="Communications: 4 pages. Full papers: no strict limit.",
    abstract_limit="200 words",
    figure_limit="Reasonable.",
    reference_format="RSC format.",
    graphical_abstract="TOC graphic.",
    submission_url="https://www.rsc.org/journals-books-databases/about-journals/dalton-transactions/",
    open_access="Hybrid.",
    apc="~£2000",
    scope="Inorganic, organometallic, bioinorganic chemistry. Catalysis, materials.",
    file_formats="Word.",
    special_notes="CIF files required for crystal structures. CCDC deposition mandatory.")

_jg("Chemistry - A European Journal",
    aliases=["chem. eur. j.", "cej", "chem eur j"],
    publisher="Wiley-VCH", issn="0947-6539",
    citation_style="Wiley (numbered).",
    word_limit="Communications: 2500 words. Full papers: no strict limit.",
    abstract_limit="200 words",
    figure_limit="Reasonable.",
    reference_format="Wiley format: numbered.",
    graphical_abstract="Required.",
    submission_url="https://chemistry-europe.onlinelibrary.wiley.com/journal/15213765",
    open_access="Hybrid.",
    apc="~€3500",
    scope="All areas of chemistry.",
    file_formats="Word.",
    special_notes="Part of Chemistry Europe family of journals.")

_jg("ACS Nano",
    aliases=["acs nano"],
    publisher="ACS", issn="1936-0851",
    citation_style="ACS style.",
    word_limit="Articles: no strict limit. Letters: ~3000 words.",
    abstract_limit="250 words",
    figure_limit="Reasonable.",
    reference_format="ACS style.",
    graphical_abstract="TOC graphic required.",
    submission_url="https://pubs.acs.org/journal/ancac3",
    open_access="Hybrid.",
    apc="~$5000",
    scope="Nanoscience and nanotechnology — synthesis, assembly, properties, applications.",
    file_formats="Word or LaTeX.",
    special_notes="Strong emphasis on characterization (TEM, AFM, DLS, etc.).")


def lookup_journal_guide(query: str) -> dict | list[dict] | None:
    """Look up journal formatting guide by name or alias."""
    q = query.strip().lower()
    if q in JOURNAL_GUIDES:
        return JOURNAL_GUIDES[q]

    matches = []
    for key, guide in JOURNAL_GUIDES.items():
        if q in key or q in guide.get("publisher", "").lower() or q in guide.get("scope", "").lower():
            if guide not in matches:
                matches.append(guide)
    return matches if matches else None


# =============================================================================
# Supporting Information (SI) checklist
# =============================================================================

SI_REQUIREMENTS: dict[str, dict] = {
    "1h_nmr": {
        "name": "¹H NMR Spectrum",
        "required_info": "Frequency (MHz), solvent, chemical shifts (δ in ppm), multiplicity (s, d, t, q, m, dd, dt, etc.), coupling constants (J in Hz), number of protons, assignment.",
        "format": "δ X.XX (mult, J = X.X Hz, NH, assignment)",
        "example": "δ 7.42 (dd, J = 8.2, 1.5 Hz, 2H, ArH), 3.85 (s, 3H, OCH₃)",
        "common_mistakes": "Missing coupling constants. Wrong multiplicity for overlapping signals. Using 'br s' without explanation. Missing residual solvent peak assignment.",
        "spectrum_required": True,
        "checklist": [
            "All peaks assigned",
            "Coupling constants for non-singlets",
            "Correct integration ratios",
            "Residual solvent peak identified",
            "Spectrum clean (no impurities >5%)",
            "Baseline flat",
        ],
    },
    "13c_nmr": {
        "name": "¹³C NMR Spectrum",
        "required_info": "Frequency (MHz), solvent, chemical shifts (δ in ppm). DEPT or HSQC data for multiplicity assignment recommended.",
        "format": "δ X.X (Cₓ assignment, optional)",
        "example": "δ 170.2 (C=O), 138.5 (C-Ar), 128.3 (CH-Ar), 55.2 (OCH₃)",
        "common_mistakes": "Missing peaks (check molecule has correct number of unique carbons). Quaternary carbons sometimes weak/missing.",
        "spectrum_required": True,
        "checklist": [
            "Number of peaks matches expected unique carbons",
            "All peaks listed in text",
            "Solvent peaks identified (CDCl₃: 77.16; DMSO-d₆: 39.52)",
            "Spectrum shows adequate S/N for all peaks",
        ],
    },
    "19f_nmr": {
        "name": "¹⁹F NMR Spectrum",
        "required_info": "Frequency (MHz), solvent, chemical shifts (δ in ppm referenced to CFCl₃ or internal standard), multiplicity, coupling constants.",
        "format": "δ -X.X (mult, J = X.X Hz)",
        "example": "δ −62.3 (s, 3F, CF₃), −110.5 (dd, J = 10.2, 8.5 Hz, 1F, ArF)",
        "spectrum_required": True,
        "checklist": ["Reference standard stated", "All F-containing groups assigned"],
    },
    "31p_nmr": {
        "name": "³¹P NMR Spectrum",
        "required_info": "Frequency (MHz), solvent, chemical shifts (δ referenced to H₃PO₄ = 0 ppm), multiplicity.",
        "format": "δ X.X (mult)",
        "example": "δ 26.5 (s)",
        "spectrum_required": True,
        "checklist": ["Reference standard stated", "Decoupled and/or coupled spectra"],
    },
    "2d_nmr": {
        "name": "2D NMR (COSY, HSQC, HMBC, NOESY)",
        "required_info": "Type of 2D experiment, frequency, solvent. Cross-peaks should be annotated.",
        "format": "HSQC, HMBC, or NOESY cross-peaks listed or annotated on spectrum.",
        "spectrum_required": True,
        "checklist": [
            "2D spectra included for structure elucidation of novel complex structures",
            "Key cross-peaks labeled",
            "Axis labels (¹H and ¹³C/¹H chemical shift axes)",
        ],
    },
    "ir": {
        "name": "Infrared Spectroscopy",
        "required_info": "Method (KBr pellet, ATR, film), key absorptions in cm⁻¹ with assignment.",
        "format": "IR (ATR): ν̃ = XXXX, XXXX, XXXX cm⁻¹",
        "example": "IR (ATR): ν̃ = 3350 (br, O-H), 1720 (s, C=O), 1600, 1510 (Ar C=C) cm⁻¹",
        "spectrum_required": True,
        "checklist": ["Key functional group absorptions listed", "Method stated"],
    },
    "hrms": {
        "name": "High-Resolution Mass Spectrometry",
        "required_info": "Ionization method (ESI, EI, APCI, MALDI), ion type ([M+H]⁺, [M+Na]⁺, [M-H]⁻), calculated mass, found mass, molecular formula.",
        "format": "HRMS (method): m/z calcd for C₁₂H₁₅NO₃ [M+ion]⁺ XXX.XXXX, found XXX.XXXX.",
        "example": "HRMS (ESI): m/z calcd for C₁₂H₁₆NO₃ [M+H]⁺ 222.1125, found 222.1128.",
        "common_mistakes": "Wrong molecular formula (forgetting to add H for [M+H]⁺). Mass accuracy >5 ppm. Wrong isotope pattern.",
        "spectrum_required": True,
        "checklist": [
            "Correct molecular formula including ion",
            "Mass accuracy within 5 ppm",
            "Ionization method stated",
            "Spectrum shows isotope pattern consistent with formula",
        ],
    },
    "melting_point": {
        "name": "Melting Point",
        "required_info": "Range (onset-endset), method (capillary, DSC), corrected/uncorrected.",
        "format": "mp X-Y °C (uncorrected)",
        "example": "mp 152-154 °C (uncorrected)",
        "spectrum_required": False,
        "checklist": ["Range not too broad (<3 °C for pure compounds)", "Method stated"],
    },
    "optical_rotation": {
        "name": "Optical Rotation",
        "required_info": "Specific rotation [α], wavelength (D-line, 589 nm), temperature, concentration (g/100 mL), solvent.",
        "format": "[α]²⁵_D = +/-X.X (c = Y.Y, solvent)",
        "example": "[α]²⁵_D = −45.2 (c = 1.0, CHCl₃)",
        "spectrum_required": False,
        "checklist": ["Temperature stated", "Concentration and solvent stated", "Sign and magnitude reported"],
    },
    "hplc": {
        "name": "HPLC Chromatogram",
        "required_info": "Column type and dimensions, mobile phase (gradient or isocratic), flow rate, detection wavelength, retention time, purity.",
        "format": "HPLC: tR = X.X min (purity: XX.X%, Column, gradient, λ = XXX nm)",
        "example": "HPLC: tR = 12.3 min (purity: 99.2%, Phenomenex Luna C18 250×4.6 mm, H₂O/MeCN 0.1% TFA, 10-90% over 30 min, 1.0 mL/min, λ = 254 nm)",
        "spectrum_required": True,
        "checklist": [
            "Full method reported (column, eluent, gradient, flow, detector)",
            "Purity stated (≥95% for biological testing)",
            "Retention time stated",
            "Chromatogram included in SI",
        ],
    },
    "chiral_hplc": {
        "name": "Chiral HPLC",
        "required_info": "Chiral column, eluent, flow rate, detection, retention times for both enantiomers, ee%.",
        "format": "Chiral HPLC: tR(major) = X.X min, tR(minor) = X.X min, ee = XX% (column, eluent, flow)",
        "spectrum_required": True,
        "checklist": [
            "Both enantiomers resolved (even if minor not detected)",
            "Racemic reference shown",
            "ee% calculated correctly",
        ],
    },
    "xray": {
        "name": "X-ray Crystallography",
        "required_info": "Crystal data, data collection parameters, structure solution/refinement details, R-factors, CCDC number.",
        "format": "CIF file deposited with CCDC (deposition number XXXXXXX).",
        "spectrum_required": False,
        "checklist": [
            "CIF file deposited with CCDC",
            "CCDC deposition number stated in manuscript",
            "ORTEP or ellipsoid plot included",
            "Crystal data table (a, b, c, α, β, γ, space group, Z, R₁, wR₂)",
            "CheckCIF alerts addressed",
            "Hydrogen atoms located or placed in calculated positions",
        ],
    },
    "elemental_analysis": {
        "name": "Elemental Analysis (CHN)",
        "required_info": "Calculated and found percentages for C, H, N (and other elements). Must agree within 0.4%.",
        "format": "Anal. Calcd for C₁₂H₁₅NO₃: C, 65.14; H, 6.83; N, 6.33. Found: C, 65.02; H, 6.79; N, 6.28.",
        "spectrum_required": False,
        "checklist": [
            "Calcd and found within 0.4% for each element",
            "Molecular formula includes any solvate/salt",
        ],
    },
    "uv_vis": {
        "name": "UV-Vis Spectroscopy",
        "required_info": "Solvent, wavelength of maxima (λmax in nm), extinction coefficients (ε in M⁻¹cm⁻¹ or L mol⁻¹ cm⁻¹).",
        "format": "UV-Vis (solvent): λmax (ε) = XXX nm (XXXX)",
        "example": "UV-Vis (MeCN): λmax (ε) = 345 nm (12,500 M⁻¹cm⁻¹), 420 nm (8,200)",
        "spectrum_required": True,
        "checklist": ["Solvent stated", "ε values calculated from Beer-Lambert", "Concentration stated"],
    },
    "fluorescence": {
        "name": "Fluorescence Spectroscopy",
        "required_info": "Excitation wavelength, emission wavelength, quantum yield (Φ), solvent, concentration.",
        "format": "Fluorescence (solvent): λex = XXX nm, λem = XXX nm, Φ = X.XX",
        "spectrum_required": True,
        "checklist": ["Reference standard for Φ stated", "Excitation and emission wavelengths"],
    },
    "tga": {
        "name": "Thermogravimetric Analysis (TGA)",
        "required_info": "Heating rate, atmosphere (N₂ or air), temperature range, mass loss events.",
        "spectrum_required": True,
        "checklist": ["Heating rate stated", "Atmosphere stated", "Onset temperatures identified"],
    },
    "dsc": {
        "name": "Differential Scanning Calorimetry (DSC)",
        "required_info": "Heating/cooling rate, temperature range, Tg, Tm, Tc values.",
        "spectrum_required": True,
        "checklist": ["Heating/cooling rates stated", "Cycle number stated (1st, 2nd heating)"],
    },
    "gc_ms": {
        "name": "GC-MS",
        "required_info": "Column, temperature program, carrier gas, injection mode, MS ionization.",
        "format": "GC-MS: tR = X.X min, m/z (% relative intensity) = XXX [M]⁺ (XX), XXX (100).",
        "spectrum_required": True,
        "checklist": ["Method details stated", "Major fragments assigned", "M⁺ peak identified"],
    },
}


def get_si_checklist(
    content_types: list[str] | None = None,
    compound_type: str = "small molecule",
) -> dict:
    """
    Generate a Supporting Information checklist.

    Args:
        content_types: List of analytical methods (e.g., ['1h_nmr', '13c_nmr', 'hrms', 'hplc']).
                       If None, returns the standard minimum for the compound type.
        compound_type: 'small molecule', 'peptide', 'polymer', 'material', 'natural product'

    Returns dict with checklist items and formatting guidance.
    """
    # Standard minimums by compound type
    standard_minimums = {
        "small molecule": ["1h_nmr", "13c_nmr", "hrms", "melting_point"],
        "peptide": ["hplc", "hrms", "1h_nmr"],
        "polymer": ["1h_nmr", "13c_nmr", "tga", "dsc"],
        "material": ["1h_nmr", "tga", "dsc", "uv_vis"],
        "natural product": ["1h_nmr", "13c_nmr", "2d_nmr", "hrms", "optical_rotation", "ir"],
    }

    if content_types is None:
        content_types = standard_minimums.get(compound_type.lower(), standard_minimums["small molecule"])

    items = []
    for ct in content_types:
        ct_lower = ct.lower().replace(" ", "_").replace("-", "_")
        if ct_lower in SI_REQUIREMENTS:
            items.append(SI_REQUIREMENTS[ct_lower])
        else:
            # Fuzzy match
            for key, val in SI_REQUIREMENTS.items():
                if ct_lower in key or ct_lower in val["name"].lower():
                    items.append(val)
                    break

    return {
        "compound_type": compound_type,
        "requested_content": content_types,
        "num_items": len(items),
        "checklist": items,
        "general_tips": [
            "Number compounds sequentially (1, 2a, 2b, 3, ...)",
            "Include General Information section (reagent sources, instrument models)",
            "Provide spectra as images (not raw data) unless journal requests FIDs",
            "Label all spectra with compound number and solvent",
            "Include a table of contents for the SI document",
            "Purity: ≥95% for compounds submitted for biological testing",
        ],
    }


# =============================================================================
# Standard chemistry abbreviations
# =============================================================================

ABBREVIATIONS: dict[str, dict[str, str]] = {
    "solvents": {
        "ACN": "acetonitrile (MeCN)",
        "DCE": "1,2-dichloroethane",
        "DCM": "dichloromethane (CH₂Cl₂)",
        "DMA": "N,N-dimethylacetamide",
        "DME": "1,2-dimethoxyethane",
        "DMF": "N,N-dimethylformamide",
        "DMSO": "dimethyl sulfoxide",
        "EtOAc": "ethyl acetate",
        "EtOH": "ethanol",
        "Et₂O": "diethyl ether",
        "Hex": "hexane(s)",
        "MeCN": "acetonitrile",
        "MeOH": "methanol",
        "NMP": "N-methyl-2-pyrrolidone",
        "PE": "petroleum ether",
        "iPrOH": "2-propanol (isopropanol)",
        "THF": "tetrahydrofuran",
        "TFE": "2,2,2-trifluoroethanol",
        "tol": "toluene",
    },
    "reagents": {
        "AIBN": "azobisisobutyronitrile",
        "BHT": "butylated hydroxytoluene",
        "Boc": "tert-butyloxycarbonyl",
        "Boc₂O": "di-tert-butyl dicarbonate",
        "BPO": "benzoyl peroxide",
        "Cbz": "benzyloxycarbonyl",
        "CSA": "camphorsulfonic acid",
        "DABCO": "1,4-diazabicyclo[2.2.2]octane",
        "DBU": "1,8-diazabicyclo[5.4.0]undec-7-ene",
        "DCC": "N,N′-dicyclohexylcarbodiimide",
        "DDQ": "2,3-dichloro-5,6-dicyano-1,4-benzoquinone",
        "DEAD": "diethyl azodicarboxylate",
        "DIAD": "diisopropyl azodicarboxylate",
        "DIBAL-H": "diisobutylaluminium hydride",
        "DIPEA": "N,N-diisopropylethylamine (Hünig's base)",
        "DMAP": "4-(dimethylamino)pyridine",
        "DMP": "Dess-Martin periodinane",
        "EDC": "1-ethyl-3-(3-dimethylaminopropyl)carbodiimide",
        "Fmoc": "9-fluorenylmethyloxycarbonyl",
        "HATU": "hexafluorophosphate azabenzotriazole tetramethyl uronium",
        "HBTU": "hexafluorophosphate benzotriazole tetramethyl uronium",
        "HOBt": "1-hydroxybenzotriazole",
        "IBX": "2-iodoxybenzoic acid",
        "LAH": "lithium aluminium hydride (LiAlH₄)",
        "LDA": "lithium diisopropylamide",
        "LiHMDS": "lithium bis(trimethylsilyl)amide",
        "mCPBA": "meta-chloroperoxybenzoic acid",
        "MOM": "methoxymethyl",
        "Ms": "methanesulfonyl (mesyl)",
        "NaHMDS": "sodium bis(trimethylsilyl)amide",
        "NBS": "N-bromosuccinimide",
        "NCS": "N-chlorosuccinimide",
        "NMO": "N-methylmorpholine N-oxide",
        "NOBIN": "2-amino-2′-hydroxy-1,1′-binaphthyl",
        "PDC": "pyridinium dichromate",
        "PCC": "pyridinium chlorochromate",
        "PMB": "para-methoxybenzyl",
        "PTSA": "para-toluenesulfonic acid",
        "SEM": "2-(trimethylsilyl)ethoxymethyl",
        "TBAF": "tetrabutylammonium fluoride",
        "TBAI": "tetrabutylammonium iodide",
        "TBDMS": "tert-butyldimethylsilyl (= TBS)",
        "TBS": "tert-butyldimethylsilyl",
        "TBDPS": "tert-butyldiphenylsilyl",
        "TCA": "trichloroacetic acid",
        "TEMPO": "(2,2,6,6-tetramethylpiperidin-1-yl)oxyl",
        "TES": "triethylsilyl",
        "Tf": "trifluoromethanesulfonyl (triflyl)",
        "TFA": "trifluoroacetic acid",
        "TfOH": "triflic acid",
        "TIPS": "triisopropylsilyl",
        "TMS": "trimethylsilyl",
        "Ts": "para-toluenesulfonyl (tosyl)",
        "Xantphos": "4,5-bis(diphenylphosphino)-9,9-dimethylxanthene",
    },
    "catalysts_ligands": {
        "BINAP": "2,2′-bis(diphenylphosphino)-1,1′-binaphthyl",
        "BrettPhos": "2-(dicyclohexylphosphino)-3,6-dimethoxy-2′,4′,6′-triisopropyl-1,1′-biphenyl",
        "cod": "1,5-cyclooctadiene",
        "Cy": "cyclohexyl",
        "dba": "dibenzylideneacetone",
        "DavePhos": "2-dicyclohexylphosphino-2′-(N,N-dimethylamino)biphenyl",
        "dppf": "1,1′-bis(diphenylphosphino)ferrocene",
        "dppe": "1,2-bis(diphenylphosphino)ethane",
        "dppp": "1,3-bis(diphenylphosphino)propane",
        "JohnPhos": "2-(di-tert-butylphosphino)biphenyl",
        "P(o-tol)₃": "tri(ortho-tolyl)phosphine",
        "PCy₃": "tricyclohexylphosphine",
        "Pd/C": "palladium on carbon",
        "Pd₂(dba)₃": "tris(dibenzylideneacetone)dipalladium(0)",
        "Pd(OAc)₂": "palladium(II) acetate",
        "Pd(PPh₃)₄": "tetrakis(triphenylphosphine)palladium(0)",
        "PPh₃": "triphenylphosphine",
        "RuPhos": "2-dicyclohexylphosphino-2′,6′-diisopropoxybiphenyl",
        "SPhos": "2-dicyclohexylphosphino-2′,6′-dimethoxybiphenyl",
        "XPhos": "2-dicyclohexylphosphino-2′,4′,6′-triisopropylbiphenyl",
        "IPr": "1,3-bis(2,6-diisopropylphenyl)imidazol-2-ylidene (NHC ligand)",
    },
    "spectroscopy": {
        "ATR": "attenuated total reflectance",
        "COSY": "correlated spectroscopy",
        "DEPT": "distortionless enhancement by polarization transfer",
        "DQF-COSY": "double quantum filtered COSY",
        "EI": "electron ionization",
        "ESI": "electrospray ionization",
        "FAB": "fast atom bombardment",
        "FID": "free induction decay",
        "HMBC": "heteronuclear multiple bond correlation",
        "HSQC": "heteronuclear single quantum coherence",
        "HRMS": "high-resolution mass spectrometry",
        "LRMS": "low-resolution mass spectrometry",
        "MALDI": "matrix-assisted laser desorption/ionization",
        "MS": "mass spectrometry",
        "NMR": "nuclear magnetic resonance",
        "NOESY": "nuclear Overhauser effect spectroscopy",
        "ROESY": "rotating frame Overhauser effect spectroscopy",
        "TOCSY": "total correlation spectroscopy",
        "TOF": "time-of-flight",
        "UV-Vis": "ultraviolet-visible spectroscopy",
    },
    "analytical": {
        "CD": "circular dichroism",
        "CE": "capillary electrophoresis",
        "DLS": "dynamic light scattering",
        "DSC": "differential scanning calorimetry",
        "FPLC": "fast protein liquid chromatography",
        "GC": "gas chromatography",
        "GPC": "gel permeation chromatography",
        "HPLC": "high-performance liquid chromatography",
        "ICP-MS": "inductively coupled plasma mass spectrometry",
        "ITC": "isothermal titration calorimetry",
        "LC-MS": "liquid chromatography-mass spectrometry",
        "MPLC": "medium-pressure liquid chromatography",
        "ORD": "optical rotatory dispersion",
        "PAGE": "polyacrylamide gel electrophoresis",
        "RP-HPLC": "reversed-phase HPLC",
        "SEC": "size exclusion chromatography",
        "SEM": "scanning electron microscopy",
        "SFC": "supercritical fluid chromatography",
        "SPR": "surface plasmon resonance",
        "TEM": "transmission electron microscopy",
        "TGA": "thermogravimetric analysis",
        "TLC": "thin-layer chromatography",
        "UPLC": "ultra-performance liquid chromatography",
        "XPS": "X-ray photoelectron spectroscopy",
        "XRD": "X-ray diffraction",
    },
    "general": {
        "aq": "aqueous",
        "br": "broad (NMR)",
        "cat.": "catalytic",
        "conc.": "concentrated",
        "d": "doublet (NMR)",
        "dd": "doublet of doublets (NMR)",
        "dt": "doublet of triplets (NMR)",
        "ee": "enantiomeric excess",
        "er": "enantiomeric ratio",
        "dr": "diastereomeric ratio",
        "equiv": "equivalent(s)",
        "m": "multiplet (NMR)",
        "M": "molar (mol/L)",
        "mp": "melting point",
        "MW": "molecular weight",
        "ppm": "parts per million",
        "q": "quartet (NMR)",
        "quant.": "quantitative",
        "Rf": "retention factor (TLC)",
        "RT": "room temperature",
        "s": "singlet (NMR)",
        "sat.": "saturated",
        "t": "triplet (NMR)",
        "tR": "retention time",
        "v/v": "volume per volume",
        "w/w": "weight per weight",
    },
    "biochemistry": {
        "ATP": "adenosine triphosphate",
        "BSA": "bovine serum albumin",
        "cDNA": "complementary DNA",
        "DMEM": "Dulbecco's modified Eagle medium",
        "DNA": "deoxyribonucleic acid",
        "DTT": "dithiothreitol",
        "EDTA": "ethylenediaminetetraacetic acid",
        "ELISA": "enzyme-linked immunosorbent assay",
        "FBS": "fetal bovine serum",
        "FRET": "Förster resonance energy transfer",
        "GSH": "glutathione (reduced)",
        "GSSG": "glutathione (oxidized)",
        "IC₅₀": "half-maximal inhibitory concentration",
        "Kd": "dissociation constant",
        "Ki": "inhibition constant",
        "Km": "Michaelis constant",
        "mRNA": "messenger RNA",
        "NAD⁺/NADH": "nicotinamide adenine dinucleotide (oxidized/reduced)",
        "PBS": "phosphate-buffered saline",
        "PCR": "polymerase chain reaction",
        "PDB": "Protein Data Bank",
        "PI": "propidium iodide",
        "RNA": "ribonucleic acid",
        "SAR": "structure-activity relationship",
        "SDS": "sodium dodecyl sulfate",
        "siRNA": "small interfering RNA",
        "TRIS": "tris(hydroxymethyl)aminomethane",
        "WT": "wild-type",
    },
}


def get_abbreviations(category: str = "all") -> dict:
    """
    Get standard chemistry abbreviations.

    Args:
        category: 'solvents', 'reagents', 'catalysts_ligands', 'spectroscopy',
                  'analytical', 'general', 'biochemistry', or 'all'

    Returns dict of abbreviation → full name mappings.
    """
    cat = category.strip().lower()
    if cat == "all":
        combined = {}
        for cat_name, abbrevs in ABBREVIATIONS.items():
            for k, v in abbrevs.items():
                combined[k] = f"{v} [{cat_name}]"
        return {
            "category": "all",
            "num_abbreviations": len(combined),
            "abbreviations": combined,
        }

    # Fuzzy match
    for key, data in ABBREVIATIONS.items():
        if cat in key or key in cat:
            return {
                "category": key,
                "num_abbreviations": len(data),
                "abbreviations": data,
            }

    return {"error": f"Unknown category: {category}. Available: {', '.join(ABBREVIATIONS.keys())}"}


def lookup_abbreviation(query: str) -> dict:
    """Look up what a specific abbreviation means."""
    q = query.strip()
    results = {}
    for cat, abbrevs in ABBREVIATIONS.items():
        for abbr, full in abbrevs.items():
            if q.lower() == abbr.lower() or q.lower() in abbr.lower():
                results[abbr] = {"meaning": full, "category": cat}
    if not results:
        # Reverse search: look in full names
        for cat, abbrevs in ABBREVIATIONS.items():
            for abbr, full in abbrevs.items():
                if q.lower() in full.lower():
                    results[abbr] = {"meaning": full, "category": cat}

    return {
        "query": q,
        "num_results": len(results),
        "results": results,
    }


# =============================================================================
# Thesis section writing guide
# =============================================================================

THESIS_GUIDES: dict[str, dict] = {
    "abstract": {
        "name": "Abstract",
        "purpose": "Concise summary of the entire thesis/paper. Should stand alone.",
        "structure": [
            "Background / Context (1-2 sentences)",
            "Research gap / Objective (1 sentence)",
            "Methods / Approach (1-2 sentences)",
            "Key results (2-3 sentences)",
            "Significance / Conclusions (1 sentence)",
        ],
        "word_limit": "Thesis: 300-500 words. Papers: 150-250 words (journal-specific).",
        "tips": [
            "Write LAST, after all other sections",
            "No references in the abstract",
            "No abbreviations without definition (or avoid altogether)",
            "Every sentence should convey essential information",
            "Include key numerical results (yields, ee's, Ki values)",
            "Use past tense for results, present tense for conclusions",
        ],
        "common_mistakes": [
            "Too vague — no specific results",
            "Too detailed — reads like a methods section",
            "Including information not in the main text",
            "Using jargon without definition",
        ],
    },
    "introduction": {
        "name": "Introduction",
        "purpose": "Establish context, identify the gap, state your contribution.",
        "structure": [
            "Broad context: What is the field? Why does it matter? (1-2 paragraphs)",
            "Literature review: What has been done? Key developments. (3-5 paragraphs)",
            "The gap: What remains unknown or unsolved? (1 paragraph)",
            "Your contribution: What did you do and why? (1 paragraph)",
            "Outline: Brief overview of thesis structure (for thesis only)",
        ],
        "tips": [
            "Funnel structure: broad → specific → your work",
            "Cite primary literature, not just reviews",
            "Be objective about others' work — acknowledge strengths and limitations",
            "Don't bury your contribution at the very end",
            "State your hypothesis clearly if applicable",
            "Use present tense for established facts, past tense for specific studies",
        ],
        "common_mistakes": [
            "Too broad introduction that reads like a textbook",
            "Literature review without clear narrative thread",
            "Gap statement is too vague ('little is known about...')",
            "Missing connection between gap and your work",
            "Citing too many reviews, not enough primary literature",
        ],
    },
    "results_discussion": {
        "name": "Results and Discussion",
        "purpose": "Present findings with interpretation. Can be combined or separate sections.",
        "structure": [
            "Organize by logical theme, not chronological order",
            "Each subsection: aim → approach → results → interpretation",
            "Figures and schemes should drive the narrative",
            "Compare with literature — agree? disagree? why?",
            "Address unexpected results honestly",
        ],
        "tips": [
            "Lead with your strongest result",
            "Every figure/table should be referenced and discussed in text",
            "Don't just describe data — interpret it",
            "Use schemes for reaction development, figures for data",
            "Discuss selectivity, scope, and limitations",
            "Compare yields/selectivities with literature benchmarks",
            "Use 'we observed' not 'it was observed' (active voice)",
        ],
        "common_mistakes": [
            "Chronological 'lab diary' organization instead of logical narrative",
            "Figures not discussed in text",
            "Over-interpreting insignificant differences",
            "Not acknowledging limitations",
            "Repeating numbers from tables in running text",
        ],
    },
    "experimental": {
        "name": "Experimental Section",
        "purpose": "Enable exact reproduction by a competent chemist.",
        "structure": [
            "General Information: instruments, reagent sources, purification of solvents",
            "General Procedures (if applicable): used for repetitive reactions",
            "Specific Compound Procedures: one per new compound",
            "Each procedure: substrate, reagent amounts (mass, mmol, equiv), conditions, workup, purification, yield",
            "Full characterization data for each new compound",
        ],
        "tips": [
            "Be specific: '5 mL' not 'some', '80 °C' not 'heated'",
            "Report actual amounts used, not theoretical",
            "Include equiv for each reagent",
            "State how reaction progress was monitored (TLC, LC-MS)",
            "Report yield as mass and percentage",
            "Use consistent formatting for all procedures",
            "Known compounds: cite literature preparation, provide ¹H NMR to confirm identity",
        ],
        "common_mistakes": [
            "Missing equivalents for reagents",
            "Vague: 'worked up in the usual manner'",
            "Inconsistent formatting between procedures",
            "Missing characterization for new compounds",
            "Not specifying which NMR solvent was used",
            "Reporting yields > 100% without explanation",
        ],
    },
    "conclusion": {
        "name": "Conclusion",
        "purpose": "Summarize key findings and their significance. Look forward.",
        "structure": [
            "Restate the problem / objective (1 sentence)",
            "Key findings — what was accomplished (2-3 paragraphs)",
            "Significance — why this matters (1 paragraph)",
            "Future work — what comes next (1 paragraph)",
        ],
        "tips": [
            "Don't just repeat the abstract",
            "Be specific about what was achieved",
            "Be honest about limitations",
            "Future work should be realistic and specific",
            "End on a forward-looking, positive note",
        ],
        "common_mistakes": [
            "Simply repeating the abstract",
            "Introducing new results not in the main text",
            "Overly speculative future work",
            "Too brief — dismissive of own work",
        ],
    },
    "supporting_information": {
        "name": "Supporting Information",
        "purpose": "Provide complete analytical data, additional experiments, and full characterization.",
        "structure": [
            "Table of Contents",
            "General Information",
            "Synthetic Procedures (if not in main text)",
            "Characterization Data (organized by compound number)",
            "NMR Spectra (¹H, ¹³C, 2D if applicable)",
            "HPLC/GC traces",
            "HRMS data",
            "X-ray crystallographic data (CIF reference)",
            "Computational details (if applicable)",
            "Additional tables and figures",
        ],
        "tips": [
            "Number compounds consistently with main text",
            "Include a table of contents for SI > 20 pages",
            "Label every spectrum with compound number and conditions",
            "Include full-page spectra, not cropped fragments",
            "Arrange spectra in compound number order",
        ],
        "common_mistakes": [
            "Unlabeled spectra",
            "Missing ¹³C NMR for new compounds",
            "HRMS without isotope pattern",
            "NMR spectra with obvious impurities not acknowledged",
        ],
    },
}


def get_thesis_guide(section: str) -> dict | None:
    """
    Get writing guidance for a thesis/paper section.

    Args:
        section: 'abstract', 'introduction', 'results_discussion', 'experimental',
                 'conclusion', 'supporting_information'
    """
    s = section.strip().lower().replace(" ", "_")
    if s in THESIS_GUIDES:
        return THESIS_GUIDES[s]

    # Fuzzy match
    for key, guide in THESIS_GUIDES.items():
        if s in key or s in guide.get("name", "").lower():
            return guide

    return {"error": f"Unknown section: {section}. Available: {', '.join(THESIS_GUIDES.keys())}"}
