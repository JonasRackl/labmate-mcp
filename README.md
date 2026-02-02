<div align="center">

# 🧪 labmate-mcp

**Your AI lab companion — from literature search to benchwork to publication.**

[![PyPI](https://img.shields.io/pypi/v/labmate-mcp?color=blue&label=PyPI)](https://pypi.org/project/labmate-mcp/)
[![Python](https://img.shields.io/badge/Python-3.10+-blue)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Tools](https://img.shields.io/badge/Tools-78-orange)](#-tool-reference)
[![Named Reactions](https://img.shields.io/badge/Named_Reactions-151-purple)](#-bench-chemistry--reference-10-tools)

**78 tools** · **25+ scientific APIs** · **151 named reactions** · **Zero config required**

[Install](#-install) · [What Can I Do?](#-what-can-i-do-with-this) · [All Tools](#-tool-reference) · [Configuration](#%EF%B8%8F-configuration) · [Examples](#-examples)

</div>

---

labmate-mcp is an [MCP server](https://modelcontextprotocol.io) that plugs into Claude and gives it deep access to scientific databases, computational chemistry, bench references, and writing utilities. **One install covers the entire research workflow.**

```
📚 Literature        ⚗️ Synthesis         🧪 Bench             📊 Analysis          ✍️ Publication
Search papers        Retrosynthesis       Named reactions      Isotope patterns     Format citations
Citation graphs      Forward prediction   Reagent calculator   Mass spectra         Build bibliography
Author profiles      Atom mapping         Protecting groups    Binding data         Experimental templates
Preprints            pKa / ADMET          Solvent reference    Crystal structures   Journal guides
Open access PDFs     NMR prediction       TLC / Column         Safety data          SI checklist
     15 tools             11 tools             27 tools             15 tools             10 tools
```

---

## 💬 What Can I Do With This?

Once installed, just talk to Claude naturally. Here are real things you can ask:

> **"Search for recent papers on copper-catalyzed C–H activation and show me the top 5 most cited"**
>
> → Searches Crossref + OpenAlex + Semantic Scholar, ranks by citations, shows abstracts

> **"I need to do a Suzuki coupling with 150 mg of my aryl bromide (MW 261). Calculate the amounts for Pd(PPh₃)₄ (5 mol%), boronic acid (1.3 eq), and K₂CO₃ (2.5 eq)"**
>
> → Returns exact masses, mmol, with the substrate as limiting reagent

> **"What's the best protecting group for a primary amine if I need it stable to acidic conditions but removable with Pd/C hydrogenation?"**
>
> → Searches 30 PGs with stability matrix, suggests Cbz

> **"Format these DOIs as an ACS-style bibliography: 10.1021/jacs.1c12345, 10.1002/anie.202112345"**
>
> → Returns numbered, journal-formatted reference list

> **"I'm writing up my experimental section. Give me a template for a Buchwald-Hartwig amination"**
>
> → Fill-in-the-blank template following journal conventions, with safety notes

> **"Look up the NMR solvent peaks for DMSO-d6 — I see something at 2.50 ppm and want to make sure it's just residual solvent"**
>
> → Residual ¹H: 2.50 ppm (quintet), ¹³C: 39.52 ppm, water: 3.33 ppm

> **"Generate a peptide library of 20 cyclic pentapeptides with at least one D-amino acid"**
>
> → Creates SMILES strings with properties (MW, logP, TPSA) for each

> **"I'm submitting to JACS. What are the formatting requirements and word limits?"**
>
> → Full guide: 5000-word limit, superscript numerals, 300-word abstract, graphical abstract specs

---

## 📦 Install

```bash
pip install labmate-mcp
```

That's it. **58 of 78 tools work immediately — no API keys, no configuration.**

### Connect to Claude

Add this to your Claude config file:

<details>
<summary><b>Claude Desktop</b> — <code>claude_desktop_config.json</code></summary>

On macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
On Windows: `%APPDATA%\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "labmate": {
      "command": "labmate-mcp"
    }
  }
}
```

</details>

<details>
<summary><b>Claude Code</b> — <code>.mcp.json</code> in your project root</summary>

```json
{
  "mcpServers": {
    "labmate": {
      "command": "labmate-mcp"
    }
  }
}
```

</details>

<details>
<summary><b>Docker</b></summary>

```bash
docker build -t labmate-mcp .
docker run -it labmate-mcp
```

</details>

Restart Claude after adding the config. You should see labmate's 78 tools available.

### Optional: Add API keys for more tools

```json
{
  "mcpServers": {
    "labmate": {
      "command": "labmate-mcp",
      "env": {
        "SEMANTIC_SCHOLAR_API_KEY": "your-key",
        "UNPAYWALL_EMAIL": "you@university.edu"
      }
    }
  }
}
```

See [Configuration](#%EF%B8%8F-configuration) for the full list.

---

## 🔧 Tool Reference

### 📚 Literature & Discovery (15 tools)

Search papers across multiple databases, explore citation graphs, find open access PDFs, and track research trends.

<details>
<summary>Show all 15 tools</summary>

| Tool | Source | What it does |
|:-----|:-------|:-------------|
| `search_papers` | Crossref + OpenAlex + S2 | Multi-source paper search with metadata fusion |
| `get_paper_details` | Crossref + OpenAlex + S2 | Full metadata: abstract, authors, citations, references |
| `find_similar_papers` | Semantic Scholar | Content-based paper recommendations |
| `get_paper_citations` | Semantic Scholar | Forward citation graph + context snippets |
| `get_paper_references` | Semantic Scholar | Backward citation graph (bibliography) |
| `get_author_profile` | OpenAlex + S2 | h-index, publications, co-authors, topics |
| `analyze_research_topic` | OpenAlex | Publication volume trends over time |
| `find_open_access_pdf` | Unpaywall | Legal open access PDF URLs |
| `search_chemrxiv` | Crossref + OpenAlex | Chemistry preprint search |
| `get_chemrxiv_categories` | — | List ChemRxiv subject categories |
| `search_web_of_science` | Web of Science | WoS search *(requires API key)* |
| `generate_bibtex` | Crossref | DOI → BibTeX (single or batch) |
| `get_journal_metrics` | OpenAlex | Impact metrics, open access %, policy |
| `search_protein_structures` | RCSB PDB | Search PDB by keyword, organism, method |
| `get_protein_structure` | RCSB PDB | Full PDB entry: resolution, ligands, sequence |

</details>

### 🔬 Compound Data & Safety (12 tools)

Look up any compound by name, SMILES, or formula. Get safety data, binding affinities, crystal structures, and more.

<details>
<summary>Show all 12 tools</summary>

| Tool | Source | What it does |
|:-----|:-------|:-------------|
| `search_compound` | PubChem | Name/SMILES/formula → compound data |
| `get_compound_properties` | PubChem | MW, SMILES, InChI, formula, XLogP, TPSA |
| `profile_compound` | Multiple | Comprehensive profile combining several databases |
| `get_safety_data` | PubChem GHS | GHS pictograms, H-statements, P-statements |
| `translate_compound_ids` | UniChem | Convert PubChem ↔ ChEMBL ↔ DrugBank ↔ ChEBI |
| `search_crystal_structures` | COD | Crystallography Open Database search |
| `search_materials_project` | Materials Project | Band gaps, formation energies *(requires key)* |
| `search_nist_webbook` | NIST | ΔHf, Cp, phase transitions, IR spectra |
| `search_mass_spectra` | MassBank | Mass spectrum search by exact mass or name |
| `search_binding_data` | BindingDB | IC₅₀, Ki, Kd binding affinities |
| `search_toxicity` | EPA CompTox | Toxicity endpoints *(requires key)* |
| `classify_natural_product` | GNPS | NP superclass / class / pathway |

</details>

### ⚗️ Computational Chemistry (11 tools)

AI-powered retrosynthesis, forward reaction prediction, pKa, solubility, ADMET, and NMR shift prediction.

<details>
<summary>Show all 11 tools</summary>

| Tool | Source | What it does |
|:-----|:-------|:-------------|
| `predict_retrosynthesis` | IBM RXN | Multi-step retrosynthetic analysis |
| `plan_synthesis` | IBM RXN | Forward synthesis route planning |
| `predict_product` | IBM RXN | Predict products from reactants + reagents |
| `predict_atom_mapping` | IBM RXN | Atom-by-atom mapping for mechanisms |
| `text_to_procedure` | IBM RXN | Natural language → structured procedure |
| `predict_pka` | Rowan Science | pKa values (any functional group, aqueous) |
| `predict_solubility` | Rowan Science | Aqueous solubility prediction |
| `predict_admet` | Rowan Science | Absorption, metabolism, toxicity prediction |
| `search_tautomers` | Rowan Science | Enumerate tautomeric forms |
| `compute_descriptors` | Rowan Science | Molecular descriptors from SMILES |
| `predict_nmr` | Rowan Science | ¹H and ¹³C chemical shift prediction |

*IBM RXN and Rowan tools require free API keys. See [Configuration](#%EF%B8%8F-configuration).*

</details>

### 🧬 Peptide Chemistry (10 tools)

Sequence-to-SMILES conversion with 450+ amino acids, cyclization, library generation, pI calculation, and MS/MS interpretation.

<details>
<summary>Show all 10 tools</summary>

| Tool | Source | What it does |
|:-----|:-------|:-------------|
| `peptide_to_smiles` | p2smi | Sequence → SMILES (450+ AAs, 5 cyclization types) |
| `peptide_cyclization_options` | p2smi | Which cyclizations does a sequence support? |
| `generate_peptide_library` | p2smi | Random peptide generation with NCAAs, D-stereo |
| `peptide_properties` | p2smi + RDKit | MW, logP, TPSA, HBD/HBA, Lipinski |
| `check_peptide_synthesis` | p2smi | SPPS feasibility: difficult motifs, aggregation |
| `modify_peptide` | p2smi | Apply N-methylation, PEGylation |
| `calculate_peptide_pi` | pichemist | Isoelectric point (8 pKa reference sets) |
| `calculate_peptide_extinction` | pep-calc.com | ε₂₈₀ (Trp/Tyr/Cys contributions) |
| `get_peptide_ion_series` | pep-calc.com | b/y/a/c/z ion ladders for MS/MS |
| `assign_peptide_ms_peaks` | pep-calc.com | Match m/z values to fragments |

*pep-calc.com tools fall back to local alternatives if the API is unavailable.*

</details>

### 🧪 Bench Chemistry — Calculators (5 tools)

The calculators every lab needs, built in. Handles unit conversions automatically.

<details>
<summary>Show all 5 tools</summary>

| Tool | What it does |
|:-----|:-------------|
| `calculate_molarity` | Solve for any unknown: mass, moles, volume, or MW |
| `calculate_dilution` | C₁V₁ = C₂V₂ with automatic unit handling |
| `calculate_reaction_mass` | Multi-reagent mass calc from equivalents |
| `calculate_yield` | Percent yield from actual / theoretical |
| `calculate_concentration` | M ↔ mM ↔ mg/mL ↔ %w/v ↔ ppm ↔ ppb |

</details>

### 📖 Bench Chemistry — Reference (10 tools)

151 named reactions. 30 protecting groups. Solvent tables. Cooling baths. TLC stains. Column guides. Buffers. NMR solvents. Everything you'd normally look up in a textbook or a dog-eared printout taped to the fume hood.

<details>
<summary>Show all 10 tools</summary>

| Tool | Coverage |
|:-----|:---------|
| `lookup_named_reaction` | **151 named reactions** with conditions, mechanism, scope, limitations |
| `lookup_protecting_group` | **30 PGs** for OH, NH, C=O, COOH with stability/lability matrix |
| `lookup_workup_procedure` | Step-by-step protocols: LAH quench, aqueous extraction, etc. |
| `lookup_solvent_properties` | **32 solvents**: bp, mp, density, polarity index, dielectric, miscibility |
| `lookup_cooling_bath` | **24 recipes** from −196 °C (lN₂) to +100 °C |
| `lookup_tlc_stain` | **13 stains** organized by functional group selectivity |
| `lookup_column_chromatography` | Solvent selection, Rf rules, loading, troubleshooting |
| `lookup_buffer_recipe` | **20+ buffers**: PBS, Tris, HEPES, TAE, TBE, RIPA, citrate, etc. |
| `lookup_amino_acid_properties` | **20 canonical AAs**: MW, pKa₁/pKa₂/pKaR, pI, hydropathy |
| `lookup_nmr_solvent` | **12 solvents**: residual ¹H/¹³C shifts, water peak, multiplicity |

</details>

### 🔧 Chemistry Utilities (5 tools)

<details>
<summary>Show all 5 tools</summary>

| Tool | What it does |
|:-----|:-------------|
| `calculate_isotope_pattern` | Isotope distribution from formula/SMILES (Cl, Br, S patterns) |
| `validate_cas_number` | CAS registry check-digit validation |
| `convert_units` | Mass, volume, energy, pressure, temperature, length, amount |
| `lookup_periodic_table` | Z, mass, electron config, electronegativity, radius, group |
| `calculate_buffer_ph` | Henderson-Hasselbalch solver with built-in pKa database |

</details>

### ✍️ Writing & Publication (10 tools)

Format citations, build bibliographies, generate experimental section templates, check journal requirements, and prepare your SI — all from within Claude.

<details>
<summary>Show all 10 tools</summary>

| Tool | Source | What it does |
|:-----|:-------|:-------------|
| `format_citation` | Crossref | DOI → formatted reference in **20+ styles** (ACS, RSC, Nature, Angew, APA…) |
| `build_bibliography` | Crossref | Batch DOIs → numbered, styled reference list |
| `lookup_iupac_name` | PubChem | SMILES → IUPAC systematic name |
| `name_to_smiles` | PubChem | Common name → SMILES + InChI + InChIKey + MW |
| `format_molecular_formula` | Local | C6H12O6 → C₆H₁₂O₆ (Unicode) / `\ce{C6H12O6}` (LaTeX) / `<sub>` (HTML) |
| `lookup_experimental_template` | Local | **18 reaction templates** with fill-in fields and safety notes |
| `lookup_journal_guide` | Local | Submission requirements for **12 top chemistry journals** |
| `generate_si_checklist` | Local | SI checklist tailored to compound type |
| `lookup_abbreviation` | Local | **193 standard abbreviations** (solvents, reagents, spectroscopy) |
| `get_thesis_guide` | Local | Section-by-section writing guide: abstract → SI |

</details>

---

## 📖 Examples

### Literature workflow

```
You:    "Find the 5 most cited papers on photoredox catalysis from 2020-2024"
Claude: [searches Crossref + Semantic Scholar, ranks by citation count]

You:    "Show me who cited paper #2 and what topics they focused on"
Claude: [retrieves forward citation graph from Semantic Scholar]

You:    "Is there an open access PDF for paper #3?"
Claude: [checks Unpaywall → returns legal OA link]

You:    "Generate BibTeX for all 5 papers"
Claude: [fetches structured citation data from Crossref]
```

### Synthesis planning

```
You:    "I want to make 4-methoxybiphenyl from 4-bromoanisole. What coupling should I use?"
Claude: [looks up Suzuki coupling conditions, suggests Pd(PPh₃)₄/K₂CO₃/PhB(OH)₂]

You:    "Calculate amounts for 200 mg scale with 5 mol% catalyst"
Claude: [reaction mass calculator → exact mg for each reagent]

You:    "What's a good workup for this?"
Claude: [aqueous workup protocol, extraction with EtOAc, Na₂SO₄ dry]
```

### Writing a paper

```
You:    "Format these 12 DOIs as an ACS bibliography"
Claude: [Crossref content negotiation → numbered reference list]

You:    "I'm writing the experimental section for a Sonogashira coupling. Give me a template"
Claude: [fill-in template with safety notes]

You:    "What SI do I need for a small molecule paper?"
Claude: [checklist: ¹H NMR, ¹³C NMR, HRMS, mp, HPLC purity, + formatting tips]

You:    "I want to submit to Angew. Chem. — what are the requirements?"
Claude: [5000-word Communication, 150-word abstract, endnote citations, graphical abstract]
```

---

## ⚙️ Configuration

**All environment variables are optional.** The server runs out of the box.

| Variable | Service | Free? | What it unlocks |
|:---------|:--------|:-----:|:----------------|
| `SEMANTIC_SCHOLAR_API_KEY` | [Semantic Scholar](https://www.semanticscholar.org/product/api#api-key) | ✅ | Higher rate limits for citations & recommendations |
| `UNPAYWALL_EMAIL` | [Unpaywall](https://unpaywall.org/) | ✅ | Open access PDF discovery |
| `RXN4CHEMISTRY_API_KEY` | [IBM RXN](https://rxn.res.ibm.com) | ✅ | Retrosynthesis, product prediction, atom mapping |
| `ROWAN_API_KEY` | [Rowan Science](https://rowan.ai) | ✅ | pKa, solubility, ADMET, tautomers, NMR prediction |
| `MATERIALS_PROJECT_API_KEY` | [Materials Project](https://materialsproject.org) | ✅ | Crystal structures, band gaps, formation energies |
| `WOS_API_KEY` | [Web of Science](https://developer.clarivate.com) | 🏛️ | Web of Science search (institutional) |
| `COMPTOX_API_KEY` | [EPA CompTox](mailto:ccte_api@epa.gov) | ✅ | Toxicity & environmental data |

**Aliases:** `S2_API_KEY` and `SEMANTIC_SCHOLAR_API_KEY` both work, as do `MP_API_KEY`/`MATERIALS_PROJECT_API_KEY` and `RXN_API_KEY`/`RXN4CHEMISTRY_API_KEY`.

---

## 🗄️ Built-in Databases

These ship with labmate and require no API calls:

| Database | Entries | What's inside |
|:---------|--------:|:--------------|
| Named reactions | **151** | Conditions, mechanism type, scope, limitations, references |
| Protecting groups | **30** | OH/NH/C=O/COOH, install/remove conditions, stability matrix |
| Solvents | **32** | bp, density, polarity index, dielectric, miscibility, safety |
| Cooling baths | **24** | Recipes from −196 °C to +100 °C |
| TLC stains | **13** | Selectivity by functional group, recipe, visualization |
| Buffer recipes | **20+** | Preparation at specific pH, temperature correction |
| Amino acids | **20** | pKa, pI, MW, hydropathy, codon, special notes |
| NMR solvents | **12** | Residual ¹H, ¹³C, water peak, multiplicity, bp |
| Experimental templates | **18** | Fill-in-the-blank for common reaction types |
| Journal guides | **12** | JACS, Angew, Nature Chem, JOC, Org Lett, etc. |
| Abbreviations | **193** | Standard abbreviations across 7 categories |
| SI requirements | **18** | Per-technique formatting and common mistakes |
| Thesis writing | **6** | Section-by-section guidance and tips |

<details>
<summary><b>All 151 named reactions</b> (click to expand)</summary>

Aldol · Appel · Arndt-Eistert · Baeyer-Villiger · Bamford-Stevens · Barton Decarboxylation · Barton-McCombie · Baylis-Hillman · Beckmann · Biginelli · Birch · Bischler-Napieralski · Bouveault-Blanc · Buchner · Buchwald-Hartwig (C–N) · Buchwald-Hartwig (C–O) · Cannizzaro · CBS · Chan-Lam · Chichibabin · Claisen Condensation · Claisen Rearrangement · Clemmensen · CuAAC Click · Comins · Cope Elimination · Corey-Bakshi-Shibata · Corey-Chaykovsky · Corey-Fuchs · Corey-Kim · Corey-Winter · Curtius · Dakin · Darzens · Dess-Martin · Dieckmann · Diels-Alder · Doering-LaFlamme · Enders SAMP/RAMP · Eschweiler-Clarke · Evans Aldol · Favorskii · Ferrier · Finkelstein · Fischer Esterification · Fischer Indole · Friedel-Crafts Acylation · Friedel-Crafts Alkylation · Fukuyama · Gabriel · Grignard · Grubbs Metathesis · Hantzsch Pyridine · Heck · Henry · Hiyama · Hiyama-Denmark · Hofmann · Horner · Horner-Wadsworth-Emmons · IBX · Ireland-Claisen · Jacobsen Epoxidation · Jones · Julia-Lythgoe · Knoevenagel · Kolbe · Kulinkovich · Kumada · Lawesson · Lemieux-Johnson · Ley · Liebeskind-Srogl · Lossen · Luche · Mannich · Matteson · Meerwein · Meerwein-Ponndorf-Verley · Michael · Midland · Minisci · Mitsunobu · Mukaiyama Aldol · Myers · Negishi · Noyori · Nozaki-Hiyama-Kishi · Ohira-Bestmann · Olefin Metathesis · Oppenauer · Oppolzer Sultam · Ozonolysis · Paal-Knorr · Parikh-Doering · Passerini · Paternò-Büchi · Petasis · Peterson · Pictet-Spengler · Pinner · Pinnick · Prins · Ramberg-Bäcklund · Reductive Amination · Reformatsky · Ring-Closing Metathesis · Ritter · Robinson Annulation · Roush · Rubottom · Sakurai-Hosomi · Sandmeyer · Sharpless AD · Sharpless AE · Shi Epoxidation · Shiina · Simmons-Smith · Sonogashira · Staudinger Ligation · Staudinger Reduction · Steglich · Stetter · Still-Gennari · Stille · Stork Enamine · Strecker · Suzuki · Suzuki-Miyaura · Swern · Takai · Tebbe · TEMPO · Tiffeneau-Demjanov · Trost AAA · Tsuji-Trost · Ugi · Ullmann · Upjohn · Vilsmeier-Haack · Wacker · Weinreb Amide · Wharton · Williamson · Wittig · Wittig Rearrangement · Wohl-Ziegler · Wolff · Wolff-Kishner · Yamaguchi · Zincke Aldehyde

</details>

---

## 🏗️ Architecture

```
labmate_mcp/
├── server.py       5,116 lines   78 MCP tool definitions + response formatting
├── bench.py        3,392 lines   Calculators + reference databases (151 reactions, 30 PGs, …)
├── apis.py         1,744 lines   HTTP clients for 25+ scientific APIs
├── writing.py      1,488 lines   Citations, templates, journal guides, SI, thesis
├── chemistry.py      572 lines   Isotope patterns, CAS, units, periodic table, pH
├── peptide.py        384 lines   p2smi + pichemist + pep-calc.com integration
└── __init__.py         4 lines   Version
                  ──────────────
                  12,700 lines
```

---

## 🤝 Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

High-impact areas: more named reactions, more experimental templates, more journal guides, tests, and bug reports.

---

## 📄 License

[MIT](LICENSE) — use freely in academia and industry.

## 📚 Citations

If labmate-mcp contributes to your research, please cite the underlying tools:

- **p2smi:** Feller, A. (2025). p2smi: Generation and analysis of drug-like peptide SMILES strings. *JOSS*, 10, 8319. [doi:10.21105/joss.08319](https://doi.org/10.21105/joss.08319)
- **pichemist:** Trastoy, B. *et al.* (2023). pIChemiSt: Structure-based isoelectric point prediction. *J. Chem. Inf. Model.* [AstraZeneca peptide-tools](https://github.com/AstraZeneca/peptide-tools)
- **OpenAlex:** Priem, J., Piwowar, H., & Orr, R. (2022). OpenAlex: A fully-open index of scholarly works. [arXiv:2205.01833](https://arxiv.org/abs/2205.01833)
