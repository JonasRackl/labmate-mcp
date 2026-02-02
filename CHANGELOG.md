# Changelog

All notable changes to labmate-mcp are documented here.

## [7.0.0] — 2025-02-02

### Added — Writing & Publication Module (10 new tools)
- `format_citation` — DOI → formatted reference in 20+ citation styles (ACS, RSC, Nature, APA, IEEE…)
- `build_bibliography` — Batch DOIs → numbered, styled reference list
- `lookup_iupac_name` — SMILES → IUPAC name via PubChem
- `name_to_smiles` — Compound name → SMILES, InChI, InChIKey via PubChem
- `format_molecular_formula` — Formula → Unicode (C₆H₁₂O₆), LaTeX, or HTML subscripts
- `lookup_experimental_template` — 18 fill-in-the-blank reaction procedure templates
- `lookup_journal_guide` — Submission formatting for 12 top chemistry journals
- `generate_si_checklist` — Supporting Information checklist by compound type
- `lookup_abbreviation` — 200+ standard chemistry abbreviations
- `get_thesis_guide` — Section-by-section thesis/paper writing guidance

### Summary
- **78 tools** spanning literature → computation → benchwork → peptides → writing
- **151 named reactions** in the reference database
- **12,700 lines** of source code across 6 modules
- Full research-to-publication workflow in a single MCP server

## [6.1.0] — 2025-02-02

### Added
- `lookup_buffer_recipe` — 20+ buffer recipes (PBS, Tris, HEPES, TAE, TBE, RIPA…)
- `lookup_amino_acid_properties` — 20 canonical amino acids with MW, pKa, pI, hydropathy
- `lookup_nmr_solvent` — 12 deuterated NMR solvents with residual shifts

## [6.0.0] — 2025-02-02

### Added — Chemistry Utilities (6 tools)
- `calculate_isotope_pattern` — Isotope distribution from formula or SMILES
- `validate_cas_number` — CAS registry check-digit validation
- `convert_units` — Scientific unit converter (mass, volume, energy, pressure, temperature)
- `lookup_periodic_table` — Element properties lookup
- `calculate_buffer_ph` — Henderson-Hasselbalch calculator
- `profile_compound` — Multi-database compound profiling

## [5.0.0] — 2025-02-02

### Added — Peptide Chemistry (10 tools)
- p2smi integration: `peptide_to_smiles`, `peptide_cyclization_options`, `generate_peptide_library`, `peptide_properties`, `check_peptide_synthesis`, `modify_peptide`
- pichemist integration: `calculate_peptide_pi`
- pep-calc.com integration: `calculate_peptide_extinction`, `get_peptide_ion_series`, `assign_peptide_ms_peaks`

### Added — Bench Chemistry (15 tools)
- 5 calculators: molarity, dilution, reaction mass, yield, concentration
- 7 reference lookups: 151 named reactions, 30 protecting groups, workup procedures, solvents, cooling baths, TLC stains, column chromatography

## [4.0.0] — 2025-02-02

### Added — Computational Chemistry (11 tools)
- IBM RXN integration: retrosynthesis, forward synthesis, product prediction, atom mapping, text-to-procedure
- Rowan Science integration: pKa, solubility, ADMET, tautomers, descriptors, NMR prediction

### Added — Compound Data (12 tools)
- PubChem, UniChem, COD, Materials Project, NIST, MassBank, BindingDB, CompTox, GHS, GNPS

## [3.0.0] — 2025-02-02

### Added — Literature Tools (15 tools)
- Crossref, OpenAlex, Semantic Scholar integration
- Citation graphs, author profiles, topic trends
- Unpaywall open access discovery
- ChemRxiv preprint search
- BibTeX generation, journal metrics
- RCSB PDB protein structure search

## [1.0.0] — 2025-02-02

Initial release with core literature search capabilities.
