# 🧪 labmate-mcp

**Your AI lab companion.** 28 tools across 19 APIs for scientific research — literature, compounds, structures, spectra, safety, synthesis, and more.

> Formerly `chemrxiv-mcp` → `scholarly-mcp` → now **labmate**, because that's what it is.

## What it does

Ask Claude anything a researcher would ask. Labmate connects to 19 free databases and returns structured answers:

- *"Search for papers on asymmetric organocatalysis"* → Semantic Scholar + OpenAlex + Crossref
- *"Get me BibTeX for these 5 DOIs"* → instant LaTeX references
- *"Tell me everything about ibuprofen"* → structure, properties, thermochemistry, safety, cross-references
- *"Find crystal structures of LiFePO4"* → COD experimental structures with CIF files
- *"Plan a retrosynthesis for aspirin"* → IBM RXN AI route planning
- *"What's the GHS hazard data for methanol?"* → pictograms, H-statements, P-statements
- *"Show me kinase inhibitor structures in PDB"* → protein-ligand complexes
- *"Where should I publish? Compare Nature Chemistry vs JACS"* → h-index, citations, OA%

## Tools (28)

### Literature & Citations (9)

| Tool | What it does |
|------|-------------|
| `search_papers` | Multi-source literature search with AI TLDRs |
| `get_paper_details` | Aggregated metadata from 4 APIs by DOI |
| `find_similar_papers` | AI-powered paper recommendations |
| `get_paper_citations` | Incoming citations with influence scores |
| `get_paper_references` | Outgoing references with influence analysis |
| `get_author_profile` | Author h-index, top papers, co-authors |
| `analyze_research_topic` | Publication trends, top authors/journals over time |
| `find_open_access_pdf` | Legal free PDF for any DOI |
| `generate_bibtex` | DOI → BibTeX (batch supported) |

### Chemistry & Compounds (10)

| Tool | What it does |
|------|-------------|
| `search_compound` | Find by name/CAS/SMILES/formula |
| `get_compound_properties` | Molecular properties & drug-likeness |
| `profile_compound` | **One query → everything** (properties + thermo + safety + cross-refs) |
| `get_safety_data` | GHS hazard pictograms, H/P statements |
| `search_nist_webbook` | Thermochemical data: ΔH, S°, Cp, bp, mp |
| `translate_compound_ids` | Cross-reference across 40+ databases |
| `classify_natural_product` | ML-powered NP pathway/class prediction |
| `search_crystal_structures` | Experimental crystal structures + CIF files |
| `search_mass_spectra` | Reference MS/MS spectra from MassBank |
| `search_binding_data` | Ki, IC50, Kd binding affinities |

### Structural Biology (2)

| Tool | What it does |
|------|-------------|
| `search_protein_structures` | Search 230k+ PDB structures |
| `get_protein_structure` | Full entry: chains, ligands, citations |

### Bibliometrics (2)

| Tool | What it does |
|------|-------------|
| `get_journal_metrics` | h-index, impact, OA%, volume trends |
| `analyze_research_topic` | Field-level publication analytics |

### Preprints (2)

| Tool | What it does |
|------|-------------|
| `search_chemrxiv` | Chemistry preprints |
| `get_chemrxiv_categories` | Subject category list |

### Optional — credential-gated (4)

| Tool | API key needed | How to get it |
|------|---------------|--------------|
| `search_materials_project` | `MP_API_KEY` | Free at [materialsproject.org](https://materialsproject.org) |
| `predict_retrosynthesis` | `RXN_API_KEY` | Free at [rxn.res.ibm.com](https://rxn.res.ibm.com) |
| `search_toxicity` | `COMPTOX_API_KEY` | Free — email ccte_api@epa.gov |
| `search_web_of_science` | `WOS_API_KEY` | Institutional subscription |

Missing keys? The tool tells you how to get one. Nothing breaks.

## Data Sources (19 APIs)

| # | API | Auth | Data |
|---|-----|------|------|
| 1 | Crossref | None | DOI metadata, BibTeX, 140M+ works |
| 2 | OpenAlex | Email | Abstracts, topics, citations, journal metrics |
| 3 | Semantic Scholar | Free key | TLDRs, recommendations, citation graphs |
| 4 | Unpaywall | Email | Open access PDF locations |
| 5 | PubChem | None | 115M+ compounds, properties, GHS safety |
| 6 | CAS Common Chemistry | None | ~500k compounds with CAS numbers |
| 7 | NIST WebBook | None* | Thermochemistry, phase data, spectra refs |
| 8 | UniChem | None | Cross-reference 40+ chemical databases |
| 9 | COD | None | 500k+ experimental crystal structures |
| 10 | MassBank EU | None | 120k+ curated MS/MS reference spectra |
| 11 | BindingDB | None | 3.2M binding measurements |
| 12 | RCSB PDB | None | 230k+ protein/nucleic acid structures |
| 13 | GNPS NPClassifier | None | Natural product ML classification |
| 14 | Materials Project | Free key | 150k+ computed materials |
| 15 | IBM RXN | Free key | AI reaction prediction & retrosynthesis |
| 16 | EPA CompTox | Free key | 1.2M chemicals with toxicity data |
| 17 | Web of Science | Institutional | Peer-reviewed literature |

*\*NIST WebBook has no official API — accessed via scraping with targeted HTML parsing.*

## Install

```bash
pip install labmate-mcp
```

### Claude Desktop config

`~/Library/Application Support/Claude/claude_desktop_config.json` (macOS):

```json
{
  "mcpServers": {
    "labmate": {
      "command": "labmate-mcp",
      "env": {
        "OPENALEX_EMAIL": "you@university.edu",
        "S2_API_KEY": "your_key",
        "UNPAYWALL_EMAIL": "you@university.edu"
      }
    }
  }
}
```

### From source

```bash
git clone https://github.com/JonasRackl/labmate-mcp.git
cd labmate-mcp
pip install -e .
```

## Configuration

All environment variables are **optional**. More keys = more features + higher rate limits.

| Variable | What it unlocks |
|----------|----------------|
| `OPENALEX_EMAIL` | 10× faster OpenAlex queries |
| `S2_API_KEY` | Dedicated Semantic Scholar rate limit |
| `UNPAYWALL_EMAIL` | Open access PDF lookups |
| `MP_API_KEY` | Materials Project tool |
| `RXN_API_KEY` | Retrosynthesis tool |
| `COMPTOX_API_KEY` | Toxicity tool |
| `WOS_API_KEY` | Web of Science tool |

## Architecture

```
labmate_mcp/
├── __init__.py      # Version
├── __main__.py      # python -m labmate_mcp
├── apis.py          # 74 API functions across 19 services (1,497 lines)
└── server.py        # 28 MCP tools with formatting (2,921 lines)
```

**Design principles:**

1. **Everything optional** — works with zero config, gets better with each key added
2. **Credential-gating** — missing API keys produce setup instructions, never errors
3. **Multi-source aggregation** — `get_paper_details` queries 4 APIs in parallel and merges results
4. **Graceful degradation** — `search_papers` falls through S2 → OpenAlex → Crossref
5. **Scraping where needed** — NIST WebBook parsed via regex (clearly documented)
6. **Orchestration** — `profile_compound` chains 6 API calls into one comprehensive report

## Migrating from scholarly-mcp / chemrxiv-mcp

```bash
pip uninstall scholarly-mcp  # or chemrxiv-mcp
pip install labmate-mcp
```

Update Claude Desktop config: change `"scholarly-mcp"` → `"labmate-mcp"` in the command field.

All existing tools work identically. 7 new tools added in v4.0.

## License

MIT
