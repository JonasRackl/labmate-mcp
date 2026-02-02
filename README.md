<div align="center">

<br>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://img.shields.io/badge/%F0%9F%A7%AA_labmate--mcp-Your_AI_Lab_Companion-0d1117?style=for-the-badge&labelColor=161b22">
  <img alt="labmate-mcp" src="https://img.shields.io/badge/%F0%9F%A7%AA_labmate--mcp-Your_AI_Lab_Companion-f6f8fa?style=for-the-badge&labelColor=eaeef2">
</picture>

### From literature search to benchwork to publication.

<br>

[![PyPI version](https://img.shields.io/pypi/v/labmate-mcp?style=flat-square&color=3572A5&label=PyPI)](https://pypi.org/project/labmate-mcp/)
[![Downloads](https://img.shields.io/pypi/dm/labmate-mcp?style=flat-square&color=3572A5&label=Downloads)](https://pypi.org/project/labmate-mcp/)
[![Python](https://img.shields.io/badge/Python-3.10+-3572A5?style=flat-square)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)](LICENSE)

<br>

**81 tools** &nbsp;·&nbsp; **25+ scientific APIs** &nbsp;·&nbsp; **202 named reactions** &nbsp;·&nbsp; **14k lines of code** &nbsp;·&nbsp; **zero config required**

<br>

[**Quick Start ↓**](#-quick-start) &nbsp;&nbsp;•&nbsp;&nbsp; [What Can I Do?](#-what-can-i-do-with-this) &nbsp;&nbsp;•&nbsp;&nbsp; [All 81 Tools](#-tool-reference) &nbsp;&nbsp;•&nbsp;&nbsp; [Configuration](#%EF%B8%8F-configuration) &nbsp;&nbsp;•&nbsp;&nbsp; [Examples](#-examples)

<br>

</div>

---

labmate-mcp is an [MCP server](https://modelcontextprotocol.io) that gives Claude deep access to scientific databases, computational chemistry, bench references, and writing utilities. **One install covers the entire research workflow** — from finding papers to writing them up.

<div align="center">

<br>

<table>
<tr>
<td align="center" width="20%">

**📚 Literature**<br>
<sub>15 tools</sub>

</td>
<td align="center" width="20%">

**⚗️ Synthesis**<br>
<sub>11 tools</sub>

</td>
<td align="center" width="20%">

**🧪 Bench**<br>
<sub>30 tools</sub>

</td>
<td align="center" width="20%">

**📊 Analysis**<br>
<sub>15 tools</sub>

</td>
<td align="center" width="20%">

**✍️ Publication**<br>
<sub>10 tools</sub>

</td>
</tr>
<tr>
<td align="center"><sub>

Search papers<br>
Citation graphs<br>
Author profiles<br>
Preprints<br>
Open access PDFs

</sub></td>
<td align="center"><sub>

Retrosynthesis<br>
Forward prediction<br>
Atom mapping<br>
pKa / ADMET<br>
NMR prediction

</sub></td>
<td align="center"><sub>

Named reactions<br>
Reagent calculator<br>
Protecting groups<br>
Solvent reference<br>
Rxn dev checklist

</sub></td>
<td align="center"><sub>

Isotope patterns<br>
Mass spectra<br>
Binding data<br>
Crystal structures<br>
Safety data

</sub></td>
<td align="center"><sub>

Format citations<br>
Build bibliography<br>
Experimental templates<br>
Journal guides<br>
SI checklist

</sub></td>
</tr>
</table>

<br>

</div>

## 🚀 Quick Start

```bash
pip install labmate-mcp
```

Then add this to your Claude config:

<details>
<summary><b>Claude Desktop</b> &nbsp;→&nbsp; <code>claude_desktop_config.json</code></summary>

<br>

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
<summary><b>Claude Code</b> &nbsp;→&nbsp; <code>.mcp.json</code> in your project root</summary>

<br>

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

<br>

```bash
docker build -t labmate-mcp .
docker run -it labmate-mcp
```

</details>

<br>

Restart Claude. **61 of 81 tools work immediately** — no API keys needed.

> [!TIP]
> Want retrosynthesis, pKa prediction, or NMR shifts? Add [free API keys](#%EF%B8%8F-configuration) for IBM RXN and Rowan Science.

---

## 💬 What Can I Do With This?

Just talk to Claude naturally. Some examples:

<table>
<tr><td>

> **"Search for recent papers on copper-catalyzed C–H activation and show me the top 5 most cited"**

Searches Crossref + OpenAlex + Semantic Scholar, ranks by citations, shows abstracts and TLDRs.

</td></tr>
<tr><td>

> **"I need to do a Suzuki coupling with 150 mg of my aryl bromide (MW 261). Calculate amounts for Pd(PPh₃)₄ (5 mol%), boronic acid (1.3 eq), and K₂CO₃ (2.5 eq)"**

Returns exact masses in mg, mmol values, and volume if a solution — with your substrate as limiting reagent.

</td></tr>
<tr><td>

> **"I'm optimising a new reaction. Walk me through the reaction development checklist"**

30-question checklist from [*Chem. Soc. Rev.* 2025](https://doi.org/10.1039/D4CS01046A): kinetics, mechanism, DoE, catalysis, scope — with 126 actionable checks and 115 practical tips.

</td></tr>
<tr><td>

> **"What's the best protecting group for a primary amine if I need it stable to acid but removable by hydrogenation?"**

Searches 30 PGs with a stability matrix (acid / base / nucleophile / oxidation / reduction / H₂-Pd) and recommends Cbz.

</td></tr>
<tr><td>

> **"Format these 12 DOIs as an ACS bibliography, then give me a Buchwald–Hartwig experimental template"**

Crossref content negotiation → numbered reference list, plus a fill-in-the-blank template with safety notes.

</td></tr>
</table>

<details>
<summary><b>More things you can ask</b></summary>

<br>

| Ask Claude… | What happens |
|:---|:---|
| "Look up the NMR solvent peaks for DMSO-d₆" | Residual ¹H: 2.50 ppm (quintet), ¹³C: 39.52 ppm, water: 3.33 ppm |
| "Generate 20 cyclic pentapeptides with D-amino acids" | SMILES strings + MW, logP, TPSA for each |
| "I'm submitting to JACS — what are the requirements?" | 5000-word limit, superscript numerals, 300-word abstract, graphical abstract specs |
| "Predict the retrosynthesis of ibuprofen" | Multi-step route from commercial starting materials |
| "What's the pKa of 4-nitrophenol?" | Quantum-chemistry prediction via Rowan Science |
| "Give me a cooling bath recipe for −42 °C" | MeCN / dry ice, or chlorobenzene / dry ice |

</details>

---

## 🔧 Tool Reference

### 📚 Literature & Discovery — 15 tools

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

### 🔬 Compound Data & Safety — 12 tools

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

### ⚗️ Computational Chemistry — 11 tools

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

### 🧬 Peptide Chemistry — 10 tools

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

</details>

### 🧪 Bench Chemistry — 18 tools

**5 calculators** for everyday lab math, plus **13 reference tools** covering 202 named reactions, a reaction development checklist, protecting groups, solvents, cooling baths, TLC stains, column guides, buffers, NMR solvents, lab tips, and safety cards.

<details>
<summary>Show all 5 calculators</summary>

| Tool | What it does |
|:-----|:-------------|
| `calculate_molarity` | Solve for any unknown: mass, moles, volume, or MW |
| `calculate_dilution` | C₁V₁ = C₂V₂ with automatic unit handling |
| `calculate_reaction_mass` | Multi-reagent mass calc from equivalents |
| `calculate_yield` | Percent yield from actual / theoretical |
| `calculate_concentration` | M ↔ mM ↔ mg/mL ↔ %w/v ↔ ppm ↔ ppb |

</details>

<details>
<summary>Show all 13 reference tools</summary>

| Tool | Coverage |
|:-----|:---------|
| `lookup_named_reaction` | **202 named reactions** — conditions, mechanism, scope, limitations |
| `lookup_rxn_dev_checklist` | **30 questions** · 126 checks · 115 tips — [Kerr *et al.*, *Chem. Soc. Rev.* 2025](https://doi.org/10.1039/D4CS01046A) |
| `lookup_protecting_group` | **30 PGs** for OH, NH, C=O, COOH with stability / lability matrix |
| `lookup_workup_procedure` | Step-by-step protocols: LAH quench, aqueous extraction, etc. |
| `lookup_solvent_properties` | **32 solvents** — bp, density, polarity index, dielectric, miscibility |
| `lookup_cooling_bath` | **24 recipes** from −196 °C (lN₂) to +100 °C |
| `lookup_tlc_stain` | **13 stains** organized by functional group selectivity |
| `lookup_column_chromatography` | Solvent selection, Rf rules, loading, troubleshooting |
| `lookup_buffer_recipe` | **20+ buffers** — PBS, Tris, HEPES, TAE, TBE, RIPA, citrate… |
| `lookup_amino_acid_properties` | **20 canonical AAs** — MW, pKa, pI, hydropathy |
| `lookup_nmr_solvent` | **12 solvents** — residual ¹H/¹³C shifts, water peak, multiplicity |
| `lookup_lab_tips` | **35 practical tips** across 9 categories |
| `lookup_safety_card` | **9 safety cards** for hazardous reagents (*n*-BuLi, NaH, LAH…) |

</details>

### 🔧 Chemistry Utilities — 5 tools

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

### ✍️ Writing & Publication — 10 tools

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

### Reaction development

```
You:    "I have a new C–H activation reaction. What should I check for mechanism?"
Claude: [rxn dev checklist → KIE, radical clocks, Hammett, Stern-Volmer, CV, DFT]

You:    "Walk me through the optimisation section"
Claude: [DoE approach, one-variable-at-a-time pitfalls, PMI/E-factor, sustainability]

You:    "What about proving catalysis vs stoichiometric?"
Claude: [mercury drop test, hot filtration, TON benchmarks, NLE analysis]
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

**Aliases:** `S2_API_KEY` / `SEMANTIC_SCHOLAR_API_KEY`, `MP_API_KEY` / `MATERIALS_PROJECT_API_KEY`, `RXN_API_KEY` / `RXN4CHEMISTRY_API_KEY` all work.

```json
{
  "mcpServers": {
    "labmate": {
      "command": "labmate-mcp",
      "env": {
        "RXN_API_KEY": "your-rxn-key",
        "ROWAN_API_KEY": "your-rowan-key",
        "UNPAYWALL_EMAIL": "you@university.edu"
      }
    }
  }
}
```

---

## 🗄️ Built-in Databases

Everything below ships with labmate — no API calls, no internet required.

<div align="center">

| | Database | Entries | What's inside |
|:--|:---------|-------:|:--------------|
| ⚗️ | Named reactions | **202** | Conditions, mechanism type, scope, limitations |
| 📋 | Rxn dev checklist | **30** questions | 126 checks + 115 tips across 7 sections |
| 🛡️ | Protecting groups | **30** | OH / NH / C=O / COOH, stability matrix |
| 🧴 | Solvents | **32** | bp, density, polarity index, dielectric, miscibility |
| ❄️ | Cooling baths | **24** | Recipes from −196 °C to +100 °C |
| 🎨 | TLC stains | **13** | Selectivity by functional group, recipe, procedure |
| 🧫 | Buffer recipes | **20+** | Preparation at specific pH, temperature correction |
| 🧬 | Amino acids | **20** | pKa, pI, MW, hydropathy, special notes |
| 📻 | NMR solvents | **12** | Residual ¹H, ¹³C, water peak, multiplicity |
| 📝 | Experimental templates | **18** | Fill-in-the-blank for common reaction types |
| 📰 | Journal guides | **12** | JACS, Angew, Nature Chem, JOC, Org Lett… |
| 🔤 | Abbreviations | **193** | Standard abbreviations across 7 categories |
| 💡 | Lab tips | **35** | Practical tips in 9 categories |
| ☣️ | Safety cards | **9** | Hazardous reagent protocols |
| 📄 | SI requirements | **18** | Per-technique formatting and common mistakes |
| 🎓 | Thesis writing | **6** | Section-by-section guidance |

</div>

<details>
<summary><b>All 202 named reactions</b></summary>

Alder-Ene · Aldol · Appel · Arbuzov · Arndt-Eistert · Baeyer-Villiger · Balz-Schiemann · Bamford-Stevens · Barton Decarboxylation · Barton-McCombie · Baylis-Hillman · Beckmann · Biginelli · Birch · Bischler-Napieralski · Blanc Chloromethylation · Bouveault-Blanc · Brown Hydroboration · Buchner Ring Expansion · Buchwald-Hartwig (C–N) · Buchwald-Hartwig (C–O) · Burgess Dehydration · Cadiot-Chodkiewicz · Cannizzaro · Carroll · Catellani · CBS · Chan-Lam · Chichibabin · Claisen Condensation · Claisen Rearrangement · Clemmensen · Click (CuAAC) · Comins · Cope Elimination · Cope Rearrangement · Corey-Bakshi-Shibata · Corey-Chaykovsky · Corey-Fuchs · Corey-Kim · Corey-Nicolaou · Corey-Winter · Cross-Metathesis · Curtius · Dakin · Darzens · Dess-Martin · Dieckmann · Diels-Alder · Doering-LaFlamme · Enders SAMP/RAMP · Eschenmoser-Claisen · Eschenmoser-Tanabe Fragmentation · Eschweiler-Clarke · Evans Aldol · Favorskii · Ferrier · Finkelstein · Fischer Esterification · Fischer Indole · Fleming-Tamao · Friedel-Crafts Acylation · Friedel-Crafts Alkylation · Fries · Fukuyama · Gabriel · Gewald · Glaser · Grignard · Grubbs Metathesis · Hantzsch Pyridine · Heck · Henry · Hiyama · Hiyama-Denmark · Hofmann · Horner · Horner-Wadsworth-Emmons · IBX · Ireland-Claisen · Jacobsen Epoxidation · Jones · Julia-Lythgoe · Kharasch · Knoevenagel · Knorr Pyrrole · Koenigs-Knorr · Kolbe · Kulinkovich · Kumada · Lawesson · Lemieux-Johnson · Ley · Liebeskind-Srogl · Lossen · Luche · Malaprade · Mander Methylenation · Mannich · Matteson · Meerwein Arylation · Meerwein Reduction · Meerwein-Ponndorf-Verley · Meinwald · Michael · Midland · Minisci · Mitsunobu · Modified Julia · Mukaiyama Aldol · Myers · Negishi · Noyori · Nozaki-Hiyama-Kishi · Ohira-Bestmann · Olefin Metathesis · Oppenauer · Oppolzer Sultam · Overman · Oxy-Cope · Ozonolysis · Paal-Knorr · Parikh-Doering · Passerini · Paternò-Büchi · Pauson-Khand · Petasis · Peterson · Pfitzner-Moffatt · Piancatelli · Pictet-Spengler · Pinner · Pinnick · Polonovski · Prevost · Prins · Ramberg-Bäcklund · Reductive Amination · Reformatsky · Rieche · Riley · Ring-Closing Metathesis · Ritter · Robinson Annulation · Roskamp · Roush · Rubottom · Saegusa-Ito · Sakurai-Hosomi · Sandmeyer · Schmidt · Shapiro · Sharpless AD · Sharpless AE · Shi Epoxidation · Shiina · Simmons-Smith · Skraup · Sonogashira · Staudinger Ligation · Staudinger Reduction · Steglich · Stetter · Still-Gennari · Stille · Stork Enamine · Strecker · Suzuki · Suzuki-Miyaura · Swern · Takai · Tebbe · TEMPO · Tiffeneau-Demjanov · Transfer Hydrogenation · Trost AAA · Tsuji-Trost · Ugi · Ullmann · Upjohn · Van Leusen · Vilsmeier-Haack · Wacker · Weinreb Amide · Wharton · Williamson · Wittig · Wittig Rearrangement · Wohl-Ziegler · Wolff · Wolff-Kishner · Yamaguchi · Zincke Aldehyde

</details>

<details>
<summary><b>Reaction development checklist — 7 sections</b></summary>

The checklist implements Kerr, Jenkinson, Sheridan & Sparr, "Reaction Development: A Student's Checklist", [*Chem. Soc. Rev.* 2025, DOI: 10.1039/D4CS01046A](https://doi.org/10.1039/D4CS01046A).

| Section | Questions | Checks | Tips |
|:--------|----------:|-------:|-----:|
| 🔍 Take Stock | 5 | 16 | 17 |
| 📈 Kinetics & Thermodynamics | 6 | 14 | 19 |
| ⚙️ Mechanism | 4 | 18 | 19 |
| 📊 Optimisation | 3 | 31 | 13 |
| 🔄 Catalysis | 4 | 16 | 16 |
| 🎯 Scope | 3 | 17 | 15 |
| 🚀 Applications | 5 | 14 | 16 |
| **Total** | **30** | **126** | **115** |

</details>

---

## 🏗️ Architecture

```
labmate_mcp/
├── server.py       5,248 lines   81 MCP tool definitions + response formatting
├── bench.py        4,714 lines   Calculators + reference databases
├── apis.py         1,744 lines   HTTP clients for 25+ scientific APIs
├── writing.py      1,488 lines   Citations, templates, journal guides, SI, thesis
├── chemistry.py      572 lines   Isotope patterns, CAS, units, periodic table, pH
├── peptide.py        384 lines   p2smi + pichemist + pep-calc.com integration
└── __init__.py         4 lines   Version
                  ──────────────
                  14,154 lines
```

---

## 🤝 Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

High-impact areas: more named reactions, more experimental templates, more journal guides, tests, and bug reports.

---

## 📄 License

[MIT](LICENSE) — use freely in academia and industry.

---

<div align="center">

## 📚 Cite

If labmate-mcp contributes to your research, please cite the underlying tools:

</div>

- **Reaction development checklist** — Kerr, M. A.; Jenkinson, M. A.; Sheridan, H.; Sparr, C. *Chem. Soc. Rev.* **2025**. [doi:10.1039/D4CS01046A](https://doi.org/10.1039/D4CS01046A)
- **p2smi** — Feller, A. *JOSS* **2025**, *10*, 8319. [doi:10.21105/joss.08319](https://doi.org/10.21105/joss.08319)
- **pichemist** — Trastoy, B. *et al.* *J. Chem. Inf. Model.* **2023**. [AstraZeneca/peptide-tools](https://github.com/AstraZeneca/peptide-tools)
- **OpenAlex** — Priem, J.; Piwowar, H.; Orr, R. [arXiv:2205.01833](https://arxiv.org/abs/2205.01833) (2022)

<div align="center">

<br>

Made with 🧪 for chemists who'd rather be in the lab than Googling.

<br>

</div>
