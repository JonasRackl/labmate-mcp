"""
labmate-mcp peptide chemistry module.

Wraps three backends:
  - p2smi (local, RDKit) — sequence → SMILES, properties, synthesis check, modifications
  - pichemist (local, AstraZeneca) — pI calculation with 8 pKa reference sets
  - pep-calc.com (REST API) — extinction coefficient, MS peak assignment, ion series

All functions are designed for MCP tool wiring: dict/str in → dict/str out.
"""

from __future__ import annotations

import json
import logging
from typing import Any

import httpx

logger = logging.getLogger(__name__)

# =============================================================================
# pep-calc.com REST API  (free, no auth, HTTP)
# =============================================================================

PEP_CALC_BASE = "http://api.pep-calc.com"


async def _pep_calc_get(endpoint: str, seq: str, n_term: str = "H", c_term: str = "OH") -> dict:
    """Call pep-calc.com API endpoint."""
    url = f"{PEP_CALC_BASE}/{endpoint}"
    params = {"seq": seq, "N_term": n_term, "C_term": c_term}
    async with httpx.AsyncClient(timeout=15) as client:
        resp = await client.get(url, params=params)
        resp.raise_for_status()
        return resp.json()


async def pep_calc_properties(seq: str, n_term: str = "H", c_term: str = "OH") -> dict:
    """Get peptide MW, formula, pI, charge summary, and extinction coefficient from pep-calc.com."""
    results = {}
    try:
        basic = await _pep_calc_get("peptide", seq, n_term, c_term)
        results["molecular_weight"] = basic.get("molecularWeight")
        results["formula"] = basic.get("formula")
        results["sequence_length"] = basic.get("seqLength")
    except Exception as e:
        results["basic_error"] = str(e)

    try:
        iso = await _pep_calc_get("peptide/iso", seq, n_term, c_term)
        results["pI"] = iso.get("pI")
    except Exception:
        pass

    try:
        charge = await _pep_calc_get("peptide/charge", seq, n_term, c_term)
        results["acidic_residues"] = charge.get("acidicCount")
        results["basic_residues"] = charge.get("basicCount")
        results["uncharged_residues"] = charge.get("unchargedCount")
    except Exception:
        pass

    try:
        ext = await _pep_calc_get("peptide/ex", seq, n_term, c_term)
        results["extinction_280nm_oxidized"] = ext.get("oxidized")
        results["extinction_280nm_reduced"] = ext.get("reduced")
    except Exception:
        pass

    return results


async def pep_calc_ms_assign(seq: str, mz_values: list[float],
                              n_term: str = "H", c_term: str = "OH") -> dict:
    """Assign observed m/z peaks to peptide deletions/modifications via pep-calc.com."""
    url = f"{PEP_CALC_BASE}/peptide/assign"
    params = {
        "seq": seq,
        "N_term": n_term,
        "C_term": c_term,
        "peaks": ",".join(str(v) for v in mz_values),
    }
    async with httpx.AsyncClient(timeout=15) as client:
        resp = await client.get(url, params=params)
        resp.raise_for_status()
        return resp.json()


async def pep_calc_ion_series(seq: str, n_term: str = "H", c_term: str = "OH") -> dict:
    """Get peptide ion series (b, y, a, c, z ions) for MS interpretation."""
    return await _pep_calc_get("peptide/ions", seq, n_term, c_term)


# =============================================================================
# pichemist (AstraZeneca) — pI calculation
# =============================================================================


def calculate_pi_from_sequence(sequence: str) -> dict:
    """
    Calculate isoelectric point using pichemist (AstraZeneca).
    Uses 8 different pKa reference sets and returns consensus pI with statistics.
    """
    from pichemist.api import pichemist_from_dict, PKaMethod
    from pichemist.model import InputAttribute

    input_dict = {
        1: {
            InputAttribute.MOL_NAME.value: sequence,
            InputAttribute.MOL_OBJECT.value: None,
            "fasta": sequence,
        }
    }
    raw = pichemist_from_dict(input_dict, method=PKaMethod.PKA_MATCHER.value)
    entry = raw.get("1", raw.get(1, {}))

    pI_data = entry.get("pI", {})
    charge_data = entry.get("QpH7", {})

    return {
        "sequence": sequence,
        "pI_mean": pI_data.get("pI mean"),
        "pI_std": pI_data.get("std"),
        "pI_stderr": pI_data.get("err"),
        "pI_interval": entry.get("pI_interval"),
        "pI_interval_threshold": entry.get("pI_interval_threshold"),
        "charge_at_pH7_mean": charge_data.get("Q at pH7.4 mean"),
        "charge_at_pH7_std": charge_data.get("std"),
        "pI_by_method": {
            k: v for k, v in pI_data.items()
            if k not in ("pI mean", "std", "err")
        },
        "charge_at_pH7_by_method": {
            k: v for k, v in charge_data.items()
            if k not in ("Q at pH7.4 mean", "std", "err")
        },
        "reference_pka_set": entry.get("pKa_set"),
    }


def calculate_pi_from_smiles(smiles: str) -> dict:
    """
    Calculate pI from a SMILES string (for modified/noncanonical peptides).
    pichemist cuts amide bonds, matches known fragments, calculates pKas for unknowns.
    """
    from pichemist.api import pichemist_from_dict, PKaMethod
    from pichemist.model import InputAttribute
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": f"Invalid SMILES: {smiles}"}

    input_dict = {
        1: {
            InputAttribute.MOL_NAME.value: smiles,
            InputAttribute.MOL_OBJECT.value: mol,
            "fasta": None,
        }
    }
    raw = pichemist_from_dict(input_dict, method=PKaMethod.PKA_MATCHER.value)
    entry = raw.get("1", raw.get(1, {}))

    pI_data = entry.get("pI", {})
    charge_data = entry.get("QpH7", {})

    return {
        "smiles": smiles,
        "pI_mean": pI_data.get("pI mean"),
        "pI_std": pI_data.get("std"),
        "pI_interval": entry.get("pI_interval"),
        "charge_at_pH7_mean": charge_data.get("Q at pH7.4 mean"),
        "pI_by_method": {
            k: v for k, v in pI_data.items()
            if k not in ("pI mean", "std", "err")
        },
    }


# =============================================================================
# p2smi — peptide SMILES generation, properties, synthesis, modifications
# =============================================================================


def sequence_to_smiles(
    sequence: str,
    cyclization: str = "",
) -> dict:
    """
    Convert peptide sequence to SMILES string.

    Args:
        sequence: Amino acid sequence (1-letter codes). Supports 450+ amino acids
                  including noncanonical (SwissSidechain).
        cyclization: One of '', 'SS', 'HT', 'SCNT', 'SCCT', 'SCSC', or
                     a manual constraint pattern like 'SSXXXCXXXCX'.
    """
    from p2smi.utilities.smilesgen import (
        linear_peptide_smiles,
        constrained_peptide_smiles,
        what_constraints,
    )

    result: dict[str, Any] = {"sequence": sequence, "cyclization": cyclization}

    if not cyclization:
        smiles = linear_peptide_smiles(sequence)
        result["smiles"] = smiles
        result["type"] = "linear"
    else:
        out = constrained_peptide_smiles(sequence, cyclization)
        if isinstance(out, tuple) and len(out) >= 3:
            result["smiles"] = out[2]
            result["applied_constraint"] = out[1]
            result["type"] = "cyclic"
        else:
            result["smiles"] = str(out)
            result["type"] = "cyclic"

    return result


def get_cyclization_options(sequence: str) -> dict:
    """Check which cyclization types a peptide sequence supports."""
    from p2smi.utilities.smilesgen import what_constraints

    raw = what_constraints(sequence)
    options = []
    for item in raw:
        if isinstance(item, (list, tuple)) and len(item) >= 2:
            constraint = item[1]
            tag = constraint[:2].upper()
            type_map = {
                "SS": "Disulfide",
                "HT": "Head-to-tail",
                "SC": "Sidechain"
            }
            # Determine subtype from mask
            mask = constraint[2:] if len(constraint) > 2 else ""
            if tag == "SC":
                has_n = "N" in mask
                has_z = "Z" in mask
                if has_n and has_z:
                    subtype = "SCSC (sidechain–sidechain)"
                elif has_n:
                    subtype = "SCNT (sidechain–N-terminus)"
                elif has_z:
                    subtype = "SCCT (sidechain–C-terminus)"
                else:
                    subtype = "SC (unspecified)"
            else:
                subtype = type_map.get(tag, tag)

            options.append({
                "type": subtype,
                "constraint_pattern": constraint,
            })

    return {
        "sequence": sequence,
        "length": len(sequence),
        "cyclization_options": options,
        "num_options": len(options),
    }


def generate_peptides(
    num_sequences: int = 10,
    min_length: int = 8,
    max_length: int = 20,
    noncanonical_percent: float = 0.0,
    dextro_percent: float = 0.0,
    cyclization: str = "none",
) -> dict:
    """
    Generate random peptide sequences with defined constraints.

    Args:
        cyclization: 'none', 'all', or comma-separated types: 'SS,HT,SCSC,SCNT,SCCT'
    """
    from p2smi.genPeps import generate_sequences

    if cyclization.lower() == "all":
        constraints = ["SS", "HT", "SCNT", "SCCT", "SCSC"]
    elif cyclization.lower() in ("none", ""):
        constraints = []
    else:
        constraints = [c.strip().upper() for c in cyclization.split(",")]

    sequences = generate_sequences(
        num_sequences=min(num_sequences, 100),
        min_length=max(min_length, 2),
        max_length=min(max_length, 100),
        noncanonical_percent=max(0.0, min(1.0, noncanonical_percent)),
        dextro_percent=max(0.0, min(1.0, dextro_percent)),
        constraints=constraints,
    )

    results = []
    for name, seq in sequences.items():
        results.append({"id": name, "sequence": seq, "length": len(seq)})

    return {
        "num_generated": len(results),
        "parameters": {
            "min_length": min_length,
            "max_length": max_length,
            "noncanonical_percent": noncanonical_percent,
            "dextro_percent": dextro_percent,
            "cyclization": cyclization,
        },
        "sequences": results,
    }


def get_peptide_properties(smiles: str) -> dict:
    """
    Compute molecular properties from a peptide SMILES string.
    Returns MW, formula, logP, TPSA, HBD, HBA, rotatable bonds, Lipinski evaluation.
    """
    from p2smi.chemProps import molecule_summary
    return molecule_summary(smiles)


def check_synthesis_feasibility(sequence: str) -> dict:
    """
    Evaluate peptide sequence for solid-phase synthesis (SPPS) feasibility.

    Checks:
    - Forbidden motifs (consecutive Pro, DG/DP, N/Q at N-terminus)
    - Cysteine content (>2 is problematic)
    - Terminal residue issues (Pro/Cys at C-terminus)
    - Glycine runs (>4 consecutive)
    - Sequence length (>50 residues)
    - Hydrophobicity (logP check)
    - Charge distribution (need charged residue every 5 positions)
    """
    from p2smi.synthRules import collect_synthesis_issues

    issues = collect_synthesis_issues(sequence)
    return {
        "sequence": sequence,
        "length": len(sequence),
        "feasible": len(issues) == 0,
        "verdict": "PASS" if not issues else "FAIL",
        "issues": issues,
        "num_issues": len(issues),
    }


def modify_peptide_smiles(
    smiles: str,
    n_methylation: bool = False,
    pegylation: bool = False,
    methylation_fraction: float = 0.3,
) -> dict:
    """
    Apply chemical modifications to a peptide SMILES string.

    Args:
        smiles: Input peptide SMILES
        n_methylation: Apply random N-methylation
        pegylation: Apply random PEGylation
        methylation_fraction: Fraction of amide sites to N-methylate (0-1)
    """
    from p2smi.chemMods import modify_sequence, is_valid_smiles

    modified, mods = modify_sequence(
        smiles,
        do_methylate=n_methylation,
        do_pegylate=pegylation,
        nmeth_residues=max(0.0, min(1.0, methylation_fraction)),
    )

    valid = is_valid_smiles(modified)

    return {
        "original_smiles": smiles,
        "modified_smiles": modified if valid else None,
        "modifications_applied": mods,
        "valid_smiles": valid,
        "error": None if valid else "Modified SMILES failed RDKit validation",
    }
