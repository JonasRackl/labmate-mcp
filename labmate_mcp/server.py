"""
labmate-mcp v7.0.0 — 78 tools from literature search to benchwork to publication.

Literature & Discovery (15): Crossref, OpenAlex, Semantic Scholar, Unpaywall, ChemRxiv, PDB
Compound Data & Safety (12): PubChem, UniChem, COD, NIST, MassBank, BindingDB, CompTox, GHS
Computational Chemistry (11): IBM RXN (retrosynthesis, prediction), Rowan (pKa, ADMET, NMR)
Peptide Chemistry (10): p2smi (450 AAs, cyclization), pichemist (pI), pep-calc.com (MS)
Bench Calculators (5): molarity, dilution, reaction mass, yield, concentration
Bench Reference (10): 151 named reactions, 30 protecting groups, solvents, TLC, column, buffers
Chemistry Utilities (5): isotope patterns, CAS validation, unit conversion, periodic table, pH
Writing & Publication (10): citations, bibliography, IUPAC, experimental templates, journal guides, SI
"""

from __future__ import annotations

import asyncio
import logging
from typing import Optional

from pydantic import BaseModel, Field
from mcp.server.fastmcp import FastMCP

from labmate_mcp.bench import (
    # Calculators
    calc_molarity,
    calc_dilution,
    calc_reaction_mass,
    calc_yield,
    calc_concentration,
    parse_formula,
    # Reference lookups
    lookup_reaction,
    lookup_pg,
    lookup_workup,
    lookup_solvent,
    lookup_cooling_bath,
    lookup_tlc_stain,
    lookup_column_guide,
    CHROM_SOLVENT_SYSTEMS,
    PROTECTING_GROUPS,
    # v6.1 additions
    lookup_buffer as _lookup_buffer,
    lookup_amino_acid as _lookup_amino_acid,
    lookup_nmr_solvent as _lookup_nmr_solvent,
    AMINO_ACIDS,
)
from labmate_mcp.chemistry import (
    validate_cas as _validate_cas,
    calculate_isotope_pattern as _calc_isotope_pattern,
    convert_units as _convert_units,
    calculate_buffer_ph as _calc_buffer_ph,
    lookup_element as _lookup_element,
    BUFFER_PKA,
    ELEMENTS,
)
from labmate_mcp.peptide import (
    sequence_to_smiles,
    get_cyclization_options,
    generate_peptides,
    get_peptide_properties,
    check_synthesis_feasibility,
    modify_peptide_smiles,
    calculate_pi_from_sequence,
    calculate_pi_from_smiles,
    pep_calc_properties,
    pep_calc_ms_assign,
    pep_calc_ion_series,
)
from labmate_mcp.writing import (
    format_citation as _format_citation,
    build_bibliography as _build_bibliography,
    iupac_from_smiles as _iupac_from_smiles,
    smiles_from_name as _smiles_from_name,
    format_molecular_formula as _format_molecular_formula,
    lookup_experimental_template as _lookup_experimental_template,
    lookup_journal_guide as _lookup_journal_guide,
    get_si_checklist as _get_si_checklist,
    lookup_abbreviation as _lookup_abbreviation,
    get_thesis_guide as _get_thesis_guide,
    get_abbreviations as _get_abbreviations,
)

from labmate_mcp.apis import (
    # Crossref
    crossref_search,
    crossref_search_chemrxiv,
    crossref_get_work,
    # OpenAlex
    openalex_search,
    openalex_get_work,
    openalex_get_author,
    openalex_search_authors,
    openalex_group_by,
    openalex_get_topic,
    openalex_get_works,
    openalex_reconstruct_abstract,
    # Semantic Scholar
    s2_search,
    s2_get_paper,
    s2_get_citations,
    s2_get_references,
    s2_get_recommendations,
    s2_search_author,
    s2_get_author,
    # Unpaywall
    unpaywall_get,
    # PubChem
    pubchem_search_by_name,
    pubchem_search_by_smiles,
    pubchem_search_by_formula,
    pubchem_get_compound,
    pubchem_get_properties,
    pubchem_get_synonyms,
    # CAS
    cas_search,
    cas_detail,
    # Web of Science
    wos_available,
    wos_search,
    # NIST WebBook (scraping)
    nist_search,
    nist_is_compound_page,
    nist_parse_search_results,
    nist_parse_compound,
    # Materials Project
    mp_available,
    mp_search,
    mp_get_material,
    # RXN4Chemistry
    rxn_available,
    rxn_predict_reaction,
    rxn_retrosynthesis,
    rxn_paragraph_to_actions,
    rxn_predict_atom_mapping,
    rxn_synthesis_plan,
    # Rowan Science
    rowan_available,
    rowan_predict_pka,
    rowan_predict_solubility,
    rowan_predict_admet,
    rowan_search_tautomers,
    rowan_compute_descriptors,
    rowan_predict_nmr,
    # UniChem
    unichem_lookup,
    # COD
    cod_search,
    # EPA CompTox
    comptox_available,
    comptox_search,
    comptox_get_details,
    comptox_get_properties,
    comptox_get_hazard,
    # MassBank
    massbank_search,
    # BindingDB
    bindingdb_by_target,
    bindingdb_by_smiles,
    # BibTeX
    crossref_get_bibtex,
    crossref_get_bibtex_batch,
    # RCSB PDB
    pdb_search,
    pdb_get_entry,
    pdb_get_entity,
    pdb_get_ligands,
    # PubChem GHS Safety
    pubchem_get_ghs,
    parse_ghs_data,
    # GNPS NPClassifier
    gnps_classify_compound,
    # OpenAlex Sources
    openalex_get_source,
    openalex_search_sources,
    # Constants
    CHEMRXIV_CATEGORIES,
    VERSION,
)

logger = logging.getLogger("scholarly-mcp")

mcp = FastMCP("labmate-mcp")


# =============================================================================
# Formatting helpers
# =============================================================================


def _authors_str(authors: list[dict], max_authors: int = 8) -> str:
    """Format author list from various API formats."""
    if not authors:
        return "N/A"
    names = []
    for a in authors[:max_authors]:
        # Semantic Scholar format
        if "name" in a:
            names.append(a["name"])
        # Crossref format
        elif "given" in a or "family" in a:
            parts = []
            if a.get("given"):
                parts.append(a["given"])
            if a.get("family"):
                parts.append(a["family"])
            name = " ".join(parts)
            if a.get("affiliation"):
                affs = [af.get("name", "") for af in a["affiliation"] if af.get("name")]
                if affs:
                    name += f" ({'; '.join(affs)})"
            names.append(name)
        else:
            names.append(str(a))
    result = ", ".join(names)
    if len(authors) > max_authors:
        result += f" ... (+{len(authors) - max_authors} more)"
    return result


def _truncate(text: str, max_len: int = 500) -> str:
    """Truncate text with ellipsis."""
    if not text or len(text) <= max_len:
        return text or ""
    return text[:max_len].rsplit(" ", 1)[0] + "..."


def _doi_link(doi: str) -> str:
    """Format DOI as clickable link."""
    if not doi:
        return ""
    return f"https://doi.org/{doi}"


def _format_s2_paper_compact(p: dict, index: int | None = None) -> str:
    """Format a Semantic Scholar paper for search results."""
    lines = []
    prefix = f"**{index}.** " if index is not None else "**"
    title = p.get("title", "Untitled")
    lines.append(f"{prefix}**{title}**" if index is None else f"{prefix}{title}")

    parts = []
    if p.get("year"):
        parts.append(str(p["year"]))
    if p.get("venue"):
        parts.append(p["venue"])
    if p.get("citationCount") is not None:
        parts.append(f"{p['citationCount']} citations")
    if p.get("influentialCitationCount"):
        parts.append(f"{p['influentialCitationCount']} influential")
    if parts:
        lines.append(f"  {' | '.join(parts)}")

    authors = p.get("authors", [])
    if authors:
        lines.append(f"  Authors: {_authors_str(authors, max_authors=5)}")

    # DOI
    ext_ids = p.get("externalIds", {})
    doi = ext_ids.get("DOI", "")
    if doi:
        lines.append(f"  DOI: {doi}")

    # TLDR (AI-generated summary)
    tldr = p.get("tldr")
    if tldr and tldr.get("text"):
        lines.append(f"  **TLDR:** {tldr['text']}")

    # Open Access
    if p.get("isOpenAccess"):
        oa_pdf = p.get("openAccessPdf", {})
        if oa_pdf and oa_pdf.get("url"):
            lines.append(f"  📄 Open Access PDF: {oa_pdf['url']}")

    return "\n".join(lines)


def _format_crossref_item_compact(item: dict, index: int | None = None) -> str:
    """Format a Crossref work for search results."""
    lines = []
    title = ""
    if item.get("title"):
        title = item["title"][0] if isinstance(item["title"], list) else item["title"]
    title = title or "Untitled"
    prefix = f"**{index}.** " if index is not None else ""
    lines.append(f"{prefix}{title}")

    parts = []
    posted = item.get("posted") or item.get("created")
    if posted and posted.get("date-parts"):
        dp = posted["date-parts"][0]
        parts.append("-".join(str(x) for x in dp if x))
    if item.get("container-title"):
        ct = item["container-title"]
        parts.append(ct[0] if isinstance(ct, list) else ct)
    if parts:
        lines.append(f"  {' | '.join(parts)}")

    authors = item.get("author", [])
    if authors:
        lines.append(f"  Authors: {_authors_str(authors, max_authors=5)}")

    doi = item.get("DOI", "")
    if doi:
        lines.append(f"  DOI: {doi}")

    abstract = item.get("abstract", "")
    if abstract:
        # Crossref abstracts sometimes have JATS XML tags
        import re
        abstract = re.sub(r"<[^>]+>", "", abstract)
        lines.append(f"  Abstract: {_truncate(abstract, 300)}")

    return "\n".join(lines)


def _format_oa_work_compact(w: dict, index: int | None = None) -> str:
    """Format an OpenAlex work for search results."""
    lines = []
    title = w.get("title", "Untitled")
    prefix = f"**{index}.** " if index is not None else ""
    lines.append(f"{prefix}{title}")

    parts = []
    if w.get("publication_year"):
        parts.append(str(w["publication_year"]))
    if w.get("cited_by_count") is not None:
        parts.append(f"{w['cited_by_count']} citations")
    loc = w.get("primary_location", {})
    if loc and loc.get("source"):
        parts.append(loc["source"].get("display_name", ""))
    if parts:
        lines.append(f"  {' | '.join(parts)}")

    authorships = w.get("authorships", [])
    if authorships:
        names = []
        for a in authorships[:5]:
            name = a.get("author", {}).get("display_name", "")
            if name:
                names.append(name)
        if names:
            lines.append(f"  Authors: {', '.join(names)}")

    doi = w.get("doi", "")
    if doi:
        lines.append(f"  DOI: {doi.replace('https://doi.org/', '')}")

    # Abstract
    inv_idx = w.get("abstract_inverted_index")
    if inv_idx:
        abstract = openalex_reconstruct_abstract(inv_idx)
        if abstract:
            lines.append(f"  Abstract: {_truncate(abstract, 300)}")

    return "\n".join(lines)


# =============================================================================
# Tool 1: search_papers — Multi-source literature search
# =============================================================================


class SearchPapersInput(BaseModel):
    """Search academic literature across multiple databases."""
    query: str = Field(
        description="Search query. Supports keywords, boolean operators (AND/OR/NOT), exact phrases in quotes."
    )
    limit: int = Field(
        default=10, ge=1, le=50,
        description="Maximum results to return (1-50)."
    )
    year: Optional[str] = Field(
        default=None,
        description="Year filter. Single year: '2024'. Range: '2020-2025'. Open: '2020-'."
    )
    open_access_only: bool = Field(
        default=False,
        description="If true, only return papers with open access PDFs."
    )
    chemistry_only: bool = Field(
        default=False,
        description="If true, filter to Chemistry field of study."
    )


@mcp.tool()
async def search_papers(params: SearchPapersInput) -> str:
    """Search academic literature. Returns papers with titles, authors, citations, and AI-generated TLDRs.

    Uses Semantic Scholar as primary source (best relevance + TLDRs),
    falls back to OpenAlex and Crossref. Supports year filtering and
    open access filtering.
    """
    # Primary: Semantic Scholar
    fos = "Chemistry" if params.chemistry_only else None
    s2_data = await s2_search(
        query=params.query,
        limit=params.limit,
        fields_of_study=fos,
        year=params.year,
        open_access_pdf=True if params.open_access_only else None,
    )

    if s2_data and s2_data.get("data"):
        papers = s2_data["data"]
        total = s2_data.get("total", len(papers))
        lines = [f"## Search Results for '{params.query}'"]
        lines.append(f"Showing {len(papers)} of {total:,} results (via Semantic Scholar)\n")
        for i, p in enumerate(papers, 1):
            lines.append(_format_s2_paper_compact(p, index=i))
            lines.append("")
        return "\n".join(lines)

    # Fallback: OpenAlex
    oa_filters = []
    if params.year:
        if "-" in params.year:
            parts = params.year.split("-")
            if parts[0]:
                oa_filters.append(f"publication_year:>{int(parts[0])-1}")
            if len(parts) > 1 and parts[1]:
                oa_filters.append(f"publication_year:<{int(parts[1])+1}")
        else:
            oa_filters.append(f"publication_year:{params.year}")
    if params.open_access_only:
        oa_filters.append("open_access.is_oa:true")

    oa_data = await openalex_search(
        query=params.query,
        per_page=params.limit,
        filters=",".join(oa_filters) if oa_filters else None,
    )

    if oa_data and oa_data.get("results"):
        works = oa_data["results"]
        total = oa_data.get("meta", {}).get("count", len(works))
        lines = [f"## Search Results for '{params.query}'"]
        lines.append(f"Showing {len(works)} of {total:,} results (via OpenAlex)\n")
        for i, w in enumerate(works, 1):
            lines.append(_format_oa_work_compact(w, index=i))
            lines.append("")
        return "\n".join(lines)

    # Fallback: Crossref
    cr_data = await crossref_search(query=params.query, rows=params.limit)
    if cr_data and cr_data.get("message", {}).get("items"):
        items = cr_data["message"]["items"]
        total = cr_data["message"].get("total-results", len(items))
        lines = [f"## Search Results for '{params.query}'"]
        lines.append(f"Showing {len(items)} of {total:,} results (via Crossref)\n")
        for i, item in enumerate(items, 1):
            lines.append(_format_crossref_item_compact(item, index=i))
            lines.append("")
        return "\n".join(lines)

    return f"No results found for '{params.query}'."


# =============================================================================
# Tool 2: get_paper_details — Full paper by DOI (aggregated)
# =============================================================================


class GetPaperInput(BaseModel):
    """Get comprehensive paper details by DOI."""
    doi: str = Field(
        description="The DOI of the paper (e.g., '10.1038/s41586-020-2649-2')."
    )


@mcp.tool()
async def get_paper_details(params: GetPaperInput) -> str:
    """Get comprehensive details for a paper by DOI. Aggregates data from
    Crossref (metadata), OpenAlex (abstract, topics, citations), Semantic Scholar
    (TLDR, influence), and Unpaywall (open access PDF).
    """
    doi = params.doi.strip()
    # Strip URL prefix if present
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.lower().startswith(prefix.lower()):
            doi = doi[len(prefix):]

    # Fetch from all sources in parallel
    cr_task = crossref_get_work(doi)
    oa_task = openalex_get_work(doi)
    s2_task = s2_get_paper(f"DOI:{doi}")
    uw_task = unpaywall_get(doi)

    results = await asyncio.gather(
        cr_task, oa_task, s2_task, uw_task, return_exceptions=True
    )
    cr_data = results[0] if not isinstance(results[0], Exception) else None
    oa_data = results[1] if not isinstance(results[1], Exception) else None
    s2_data = results[2] if not isinstance(results[2], Exception) else None
    uw_data = results[3] if not isinstance(results[3], Exception) else None

    if not any([cr_data, oa_data, s2_data]):
        return f"Paper not found for DOI: {doi}"

    lines = []

    # --- Title ---
    title = (
        (s2_data or {}).get("title")
        or (oa_data or {}).get("title")
        or ((cr_data or {}).get("title", [""])[0] if isinstance((cr_data or {}).get("title"), list) else (cr_data or {}).get("title", ""))
        or "Untitled"
    )
    lines.append(f"# {title}")
    lines.append(f"**DOI:** {doi}")

    # --- Date ---
    date = (
        (s2_data or {}).get("publicationDate")
        or str((s2_data or {}).get("year", ""))
    )
    if not date and cr_data:
        posted = cr_data.get("posted") or cr_data.get("created") or {}
        if posted.get("date-parts"):
            dp = posted["date-parts"][0]
            date = "-".join(str(x) for x in dp if x)
    if not date and oa_data:
        date = oa_data.get("publication_date", str(oa_data.get("publication_year", "")))
    if date:
        lines.append(f"**Published:** {date}")

    # --- Venue / Journal ---
    venue = (s2_data or {}).get("venue") or ""
    if not venue and oa_data:
        loc = oa_data.get("primary_location", {})
        if loc and loc.get("source"):
            venue = loc["source"].get("display_name", "")
    if not venue and cr_data:
        ct = cr_data.get("container-title", [])
        venue = ct[0] if ct else ""
    if venue:
        lines.append(f"**Journal:** {venue}")

    # --- License ---
    license_url = ""
    if cr_data and cr_data.get("license"):
        license_url = cr_data["license"][0].get("URL", "")
    if not license_url and oa_data:
        loc = oa_data.get("best_oa_location") or oa_data.get("primary_location", {})
        if loc:
            license_url = loc.get("license", "") or ""
    if license_url:
        lines.append(f"**License:** {license_url}")

    # --- Topics ---
    topics = []
    if oa_data and oa_data.get("topics"):
        topics = [t.get("display_name", "") for t in oa_data["topics"][:5] if t.get("display_name")]
    elif s2_data and s2_data.get("s2FieldsOfStudy"):
        topics = [f.get("category", "") for f in s2_data["s2FieldsOfStudy"][:5] if f.get("category")]
    if topics:
        lines.append(f"**Topics:** {', '.join(topics)}")

    lines.append("")

    # --- Authors ---
    lines.append("## Authors")
    if cr_data and cr_data.get("author"):
        for a in cr_data["author"]:
            parts = []
            if a.get("given"):
                parts.append(a["given"])
            if a.get("family"):
                parts.append(a["family"])
            name = " ".join(parts)
            affs = [af.get("name", "") for af in a.get("affiliation", []) if af.get("name")]
            if affs:
                name += f" ({'; '.join(affs)})"
            if a.get("ORCID"):
                name += f" [ORCID]({a['ORCID']})"
            lines.append(f"  - {name}")
    elif s2_data and s2_data.get("authors"):
        for a in s2_data["authors"]:
            name = a.get("name", "Unknown")
            if a.get("authorId"):
                name += f" (S2:{a['authorId']})"
            lines.append(f"  - {name}")
    elif oa_data and oa_data.get("authorships"):
        for a in oa_data["authorships"]:
            name = a.get("author", {}).get("display_name", "Unknown")
            insts = [i.get("display_name", "") for i in a.get("institutions", []) if i.get("display_name")]
            if insts:
                name += f" ({'; '.join(insts)})"
            lines.append(f"  - {name}")

    lines.append("")

    # --- Abstract ---
    abstract = ""
    if s2_data and s2_data.get("abstract"):
        abstract = s2_data["abstract"]
    elif oa_data and oa_data.get("abstract_inverted_index"):
        abstract = openalex_reconstruct_abstract(oa_data["abstract_inverted_index"])
    elif cr_data and cr_data.get("abstract"):
        import re
        abstract = re.sub(r"<[^>]+>", "", cr_data["abstract"])
    if abstract:
        lines.append("## Abstract")
        lines.append(abstract)
        lines.append("")

    # --- TLDR (AI Summary) ---
    if s2_data and s2_data.get("tldr") and s2_data["tldr"].get("text"):
        lines.append(f"## AI Summary (TLDR)")
        lines.append(s2_data["tldr"]["text"])
        lines.append("")

    # --- Citation metrics ---
    lines.append("## Metrics")
    cite_count = (
        (s2_data or {}).get("citationCount")
        or (oa_data or {}).get("cited_by_count")
    )
    ref_count = (s2_data or {}).get("referenceCount")
    influential = (s2_data or {}).get("influentialCitationCount")

    metrics_parts = []
    if cite_count is not None:
        metrics_parts.append(f"**Citations:** {cite_count}")
    if influential:
        metrics_parts.append(f"**Influential citations:** {influential}")
    if ref_count is not None:
        metrics_parts.append(f"**References:** {ref_count}")

    # OpenAlex normalized metrics
    if oa_data:
        fwci = oa_data.get("fwci")
        if fwci is not None:
            metrics_parts.append(f"**FWCI:** {fwci:.2f}")
        percentile = oa_data.get("citation_normalized_percentile", {})
        if percentile and percentile.get("value") is not None:
            metrics_parts.append(f"**Citation percentile:** {percentile['value']:.1%}")

    if metrics_parts:
        lines.append("  |  ".join(metrics_parts))
    else:
        lines.append("No citation data available yet.")
    lines.append("")

    # --- Open Access ---
    oa_url = None
    oa_status = None

    if uw_data and not isinstance(uw_data, Exception):
        oa_status = uw_data.get("oa_status", "unknown")
        best_loc = uw_data.get("best_oa_location")
        if best_loc:
            oa_url = best_loc.get("url_for_pdf") or best_loc.get("url")

    if not oa_url and s2_data and s2_data.get("openAccessPdf"):
        oa_url = s2_data["openAccessPdf"].get("url")

    if not oa_url and oa_data:
        oa_locs = oa_data.get("open_access", {})
        if oa_locs.get("oa_url"):
            oa_url = oa_locs["oa_url"]
            oa_status = oa_locs.get("oa_status")

    if oa_url or oa_status:
        lines.append("## Open Access")
        if oa_status:
            lines.append(f"**Status:** {oa_status}")
        if oa_url:
            lines.append(f"**PDF:** {oa_url}")
        lines.append("")

    # --- Links ---
    lines.append("## Links")
    lines.append(f"- DOI: {_doi_link(doi)}")
    if s2_data and s2_data.get("url"):
        lines.append(f"- Semantic Scholar: {s2_data['url']}")
    if oa_data and oa_data.get("id"):
        lines.append(f"- OpenAlex: {oa_data['id']}")

    return "\n".join(lines)


# =============================================================================
# Tool 3: find_similar_papers — Semantic Scholar recommendations
# =============================================================================


class FindSimilarInput(BaseModel):
    """Find papers similar to a given paper."""
    doi: str = Field(
        description="DOI of the paper to find similar papers for."
    )
    limit: int = Field(
        default=10, ge=1, le=50,
        description="Number of recommendations (1-50)."
    )


@mcp.tool()
async def find_similar_papers(params: FindSimilarInput) -> str:
    """Find papers similar to a given paper using Semantic Scholar's
    recommendation engine. Based on citation patterns and content similarity.
    Returns papers with TLDRs and citation counts.
    """
    doi = params.doi.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.lower().startswith(prefix.lower()):
            doi = doi[len(prefix):]

    data = await s2_get_recommendations(
        positive_ids=[f"DOI:{doi}"],
        limit=params.limit,
    )

    if not data or not data.get("recommendedPapers"):
        return f"No recommendations found for DOI:{doi}. The paper may not be indexed in Semantic Scholar."

    papers = data["recommendedPapers"]
    lines = [f"## Papers Similar to DOI:{doi}"]
    lines.append(f"Found {len(papers)} recommendations\n")
    for i, p in enumerate(papers, 1):
        lines.append(_format_s2_paper_compact(p, index=i))
        lines.append("")
    return "\n".join(lines)


# =============================================================================
# Tool 4: get_paper_citations — Citation graph (incoming)
# =============================================================================


class GetCitationsInput(BaseModel):
    """Get papers that cite a given paper."""
    doi: str = Field(description="DOI of the paper.")
    limit: int = Field(default=20, ge=1, le=100, description="Max citations to return.")
    influential_only: bool = Field(
        default=False,
        description="If true, only return influential citations (key references, not just mentions)."
    )


@mcp.tool()
async def get_paper_citations(params: GetCitationsInput) -> str:
    """Get papers that cite a given work. Includes influence scores and
    citation intents (methodology, background, result comparison).
    Uses Semantic Scholar's citation graph.
    """
    doi = params.doi.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.lower().startswith(prefix.lower()):
            doi = doi[len(prefix):]

    data = await s2_get_citations(f"DOI:{doi}", limit=params.limit)
    if not data or not data.get("data"):
        return f"No citation data found for DOI:{doi}."

    citations = data["data"]
    if params.influential_only:
        citations = [c for c in citations if c.get("isInfluential")]

    lines = [f"## Papers Citing DOI:{doi}"]
    lines.append(f"Showing {len(citations)} citations")
    if params.influential_only:
        lines.append("*(filtered to influential citations only)*")
    lines.append("")

    for i, c in enumerate(citations, 1):
        cp = c.get("citingPaper", {})
        title = cp.get("title", "Untitled")
        year = cp.get("year", "")
        cite_count = cp.get("citationCount", 0)
        authors = _authors_str(cp.get("authors", []), max_authors=3)

        influential = "⭐ " if c.get("isInfluential") else ""
        intents = c.get("intents", [])
        intent_str = f" [{', '.join(intents)}]" if intents else ""

        lines.append(f"**{i}.** {influential}{title}")
        lines.append(f"  {year} | {cite_count} citations | {authors}{intent_str}")

        # Citation context (how the paper is cited)
        contexts = c.get("contexts", [])
        if contexts:
            lines.append(f"  Context: \"{_truncate(contexts[0], 200)}\"")
        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 5: get_paper_references — Citation graph (outgoing)
# =============================================================================


class GetReferencesInput(BaseModel):
    """Get papers referenced by a given paper."""
    doi: str = Field(description="DOI of the paper.")
    limit: int = Field(default=20, ge=1, le=100, description="Max references to return.")
    influential_only: bool = Field(
        default=False,
        description="If true, only return influential references."
    )


@mcp.tool()
async def get_paper_references(params: GetReferencesInput) -> str:
    """Get the reference list of a paper with influence analysis.
    Shows which references are influential (key vs. background citations).
    Uses Semantic Scholar's citation graph.
    """
    doi = params.doi.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.lower().startswith(prefix.lower()):
            doi = doi[len(prefix):]

    data = await s2_get_references(f"DOI:{doi}", limit=params.limit)
    if not data or not data.get("data"):
        return f"No reference data found for DOI:{doi}."

    refs = data["data"]
    if params.influential_only:
        refs = [r for r in refs if r.get("isInfluential")]

    lines = [f"## References from DOI:{doi}"]
    lines.append(f"Showing {len(refs)} references")
    if params.influential_only:
        lines.append("*(filtered to influential references only)*")
    lines.append("")

    for i, r in enumerate(refs, 1):
        rp = r.get("citedPaper", {})
        title = rp.get("title", "Untitled")
        year = rp.get("year", "")
        cite_count = rp.get("citationCount", 0)
        authors = _authors_str(rp.get("authors", []), max_authors=3)
        influential = "⭐ " if r.get("isInfluential") else ""
        intents = r.get("intents", [])
        intent_str = f" [{', '.join(intents)}]" if intents else ""

        lines.append(f"**{i}.** {influential}{title}")
        lines.append(f"  {year} | {cite_count} citations | {authors}{intent_str}")
        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 6: get_author_profile — Author metrics & publications
# =============================================================================


class GetAuthorInput(BaseModel):
    """Look up an academic author profile."""
    name: str = Field(
        description="Author name to search for (e.g., 'John Smith', 'K. Jorner')."
    )


@mcp.tool()
async def get_author_profile(params: GetAuthorInput) -> str:
    """Look up an author's profile including h-index, citation count,
    publication history, and institutional affiliation. Searches both
    Semantic Scholar and OpenAlex for comprehensive coverage.
    """
    name = params.name.strip()

    # Search both sources in parallel
    s2_task = s2_search_author(name, limit=3)
    oa_task = openalex_search_authors(name, per_page=3)

    s2_results, oa_results = await asyncio.gather(
        s2_task, oa_task, return_exceptions=True
    )

    lines = [f"## Author Profile: {name}\n"]

    # Semantic Scholar results
    if s2_results and not isinstance(s2_results, Exception) and s2_results.get("data"):
        for author in s2_results["data"][:2]:
            lines.append(f"### {author.get('name', name)} (Semantic Scholar)")
            parts = []
            if author.get("hIndex") is not None:
                parts.append(f"h-index: {author['hIndex']}")
            if author.get("citationCount") is not None:
                parts.append(f"Citations: {author['citationCount']:,}")
            if author.get("paperCount") is not None:
                parts.append(f"Papers: {author['paperCount']}")
            if parts:
                lines.append("  " + "  |  ".join(parts))
            affs = author.get("affiliations", [])
            if affs:
                lines.append(f"  Affiliations: {', '.join(affs)}")
            lines.append(f"  S2 ID: {author.get('authorId', 'N/A')}")

            # Fetch detailed papers for top result
            if author.get("authorId"):
                detail = await s2_get_author(author["authorId"])
                if detail and detail.get("papers"):
                    papers = sorted(
                        detail["papers"],
                        key=lambda p: p.get("citationCount", 0),
                        reverse=True,
                    )[:5]
                    lines.append("  **Top papers:**")
                    for p in papers:
                        title = p.get("title", "Untitled")
                        year = p.get("year", "")
                        cc = p.get("citationCount", 0)
                        doi = p.get("externalIds", {}).get("DOI", "")
                        doi_str = f" (DOI: {doi})" if doi else ""
                        lines.append(f"    - {title} ({year}, {cc} citations){doi_str}")
            lines.append("")

    # OpenAlex results
    if oa_results and not isinstance(oa_results, Exception) and oa_results.get("results"):
        for author in oa_results["results"][:2]:
            lines.append(f"### {author.get('display_name', name)} (OpenAlex)")
            parts = []
            if author.get("summary_stats"):
                stats = author["summary_stats"]
                if stats.get("h_index") is not None:
                    parts.append(f"h-index: {stats['h_index']}")
                if stats.get("i10_index") is not None:
                    parts.append(f"i10-index: {stats['i10_index']}")
            if author.get("cited_by_count") is not None:
                parts.append(f"Citations: {author['cited_by_count']:,}")
            if author.get("works_count") is not None:
                parts.append(f"Works: {author['works_count']}")
            if parts:
                lines.append("  " + "  |  ".join(parts))

            # Affiliations
            affs = author.get("affiliations", [])
            if affs:
                inst_names = [a.get("institution", {}).get("display_name", "") for a in affs[:3]]
                lines.append(f"  Affiliations: {', '.join(n for n in inst_names if n)}")

            # Last known institution
            last_inst = author.get("last_known_institutions", [])
            if last_inst:
                inst = last_inst[0]
                lines.append(f"  Current: {inst.get('display_name', 'N/A')} ({inst.get('country_code', '')})")

            if author.get("orcid"):
                lines.append(f"  ORCID: {author['orcid']}")
            lines.append(f"  OpenAlex: {author.get('id', 'N/A')}")
            lines.append("")

    if len(lines) <= 2:
        return f"No author profiles found for '{name}'."

    return "\n".join(lines)


# =============================================================================
# Tool 7: analyze_research_topic — Bibliometric trends
# =============================================================================


class AnalyzeTopicInput(BaseModel):
    """Analyze research trends for a topic."""
    topic: str = Field(
        description="Research topic to analyze (e.g., 'asymmetric catalysis', 'CRISPR', 'perovskite solar cells')."
    )
    year_from: int = Field(default=2015, description="Start year for trend analysis.")
    year_to: int = Field(default=2025, description="End year for trend analysis.")


@mcp.tool()
async def analyze_research_topic(params: AnalyzeTopicInput) -> str:
    """Analyze publication trends for a research topic. Returns:
    - Publication volume per year
    - Top authors in the field
    - Top journals publishing on this topic
    - Open access proportion

    Uses OpenAlex bibliometric aggregations.
    """
    topic = params.topic.strip()
    year_filter = f"publication_year:{params.year_from}-{params.year_to}"

    # Run multiple aggregations in parallel
    trend_task = openalex_group_by(
        group_by="publication_year",
        search=topic,
        filters=year_filter,
    )
    author_task = openalex_group_by(
        group_by="authorships.author.id",
        search=topic,
        filters=year_filter,
        per_page=15,
    )
    journal_task = openalex_group_by(
        group_by="primary_location.source.id",
        search=topic,
        filters=year_filter,
        per_page=15,
    )
    oa_task = openalex_group_by(
        group_by="open_access.oa_status",
        search=topic,
        filters=year_filter,
    )

    trend_data, author_data, journal_data, oa_data = await asyncio.gather(
        trend_task, author_task, journal_task, oa_task, return_exceptions=True
    )

    lines = [f"## Research Landscape: {topic}"]
    lines.append(f"Analysis period: {params.year_from}–{params.year_to}\n")

    # Total count
    total = 0
    if trend_data and not isinstance(trend_data, Exception):
        meta = trend_data.get("meta", {})
        total = meta.get("count", 0)
    lines.append(f"**Total publications:** {total:,}\n")

    # --- Publication trend ---
    if trend_data and not isinstance(trend_data, Exception) and trend_data.get("group_by"):
        lines.append("### Publication Trend (per year)")
        groups = sorted(trend_data["group_by"], key=lambda g: g.get("key", ""))
        for g in groups:
            year = g.get("key", "")
            count = g.get("count", 0)
            bar = "█" * max(1, int(count / max(g2.get("count", 1) for g2 in groups) * 30))
            lines.append(f"  {year}: {bar} {count:,}")
        lines.append("")

    # --- Top authors ---
    if author_data and not isinstance(author_data, Exception) and author_data.get("group_by"):
        lines.append("### Top Authors")
        for g in author_data["group_by"][:10]:
            name = g.get("key_display_name", g.get("key", "Unknown"))
            count = g.get("count", 0)
            lines.append(f"  - {name}: {count} publications")
        lines.append("")

    # --- Top journals ---
    if journal_data and not isinstance(journal_data, Exception) and journal_data.get("group_by"):
        lines.append("### Top Journals")
        for g in journal_data["group_by"][:10]:
            name = g.get("key_display_name", g.get("key", "Unknown"))
            count = g.get("count", 0)
            lines.append(f"  - {name}: {count} publications")
        lines.append("")

    # --- Open access ---
    if oa_data and not isinstance(oa_data, Exception) and oa_data.get("group_by"):
        lines.append("### Open Access Distribution")
        for g in oa_data["group_by"]:
            status = g.get("key", "unknown")
            count = g.get("count", 0)
            pct = (count / total * 100) if total > 0 else 0
            lines.append(f"  - {status}: {count:,} ({pct:.1f}%)")
        lines.append("")

    if len(lines) <= 4:
        return f"Insufficient data for topic analysis on '{topic}'."

    return "\n".join(lines)


# =============================================================================
# Tool 8: find_open_access_pdf — Unpaywall OA discovery
# =============================================================================


class FindOAInput(BaseModel):
    """Find open access PDF for a paper."""
    doi: str = Field(description="DOI of the paper.")


@mcp.tool()
async def find_open_access_pdf(params: FindOAInput) -> str:
    """Find a free, legal PDF for a paper using Unpaywall.
    Returns direct PDF link, OA status (gold/green/hybrid/bronze),
    and the repository or publisher hosting it.
    """
    doi = params.doi.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if doi.lower().startswith(prefix.lower()):
            doi = doi[len(prefix):]

    data = await unpaywall_get(doi)
    if not data:
        return f"Could not look up DOI: {doi}. Verify the DOI is correct."

    lines = [f"## Open Access Status for {doi}\n"]
    lines.append(f"**Title:** {data.get('title', 'N/A')}")
    lines.append(f"**Journal:** {data.get('journal_name', 'N/A')}")
    lines.append(f"**OA Status:** {data.get('oa_status', 'unknown')}")
    lines.append(f"**Is OA:** {'Yes' if data.get('is_oa') else 'No'}")
    lines.append("")

    if data.get("is_oa"):
        best = data.get("best_oa_location", {})
        if best:
            lines.append("### Best Open Access Location")
            if best.get("url_for_pdf"):
                lines.append(f"**📄 PDF:** {best['url_for_pdf']}")
            if best.get("url"):
                lines.append(f"**🔗 Landing page:** {best['url']}")
            if best.get("license"):
                lines.append(f"**License:** {best['license']}")
            if best.get("host_type"):
                lines.append(f"**Host:** {best['host_type']}")
            if best.get("version"):
                lines.append(f"**Version:** {best['version']}")
            lines.append("")

        # Additional locations
        all_locs = data.get("oa_locations", [])
        if len(all_locs) > 1:
            lines.append("### Other Locations")
            for loc in all_locs[1:4]:  # Show up to 3 more
                url = loc.get("url_for_pdf") or loc.get("url", "")
                host = loc.get("host_type", "")
                version = loc.get("version", "")
                lines.append(f"  - {url} ({host}, {version})")
            lines.append("")
    else:
        lines.append("*This paper is not currently available as open access.*")
        lines.append(f"Publisher page: {_doi_link(doi)}")

    return "\n".join(lines)


# =============================================================================
# Tool 9: search_chemrxiv — ChemRxiv preprints
# =============================================================================


class SearchChemRxivInput(BaseModel):
    """Search ChemRxiv chemistry preprints."""
    query: str = Field(
        default="",
        description="Search query (keywords). Leave empty for recent preprints."
    )
    limit: int = Field(default=10, ge=1, le=50, description="Max results (1-50).")
    date_from: Optional[str] = Field(
        default=None,
        description="Start date YYYY-MM-DD (e.g., '2024-01-01')."
    )
    date_to: Optional[str] = Field(
        default=None,
        description="End date YYYY-MM-DD."
    )
    sort: str = Field(
        default="relevance",
        description="Sort order: 'relevance' or 'date'."
    )


@mcp.tool()
async def search_chemrxiv(params: SearchChemRxivInput) -> str:
    """Search ChemRxiv chemistry preprints via Crossref.
    Enriches results with abstracts from OpenAlex.

    IMPORTANT: ChemRxiv preprints are NOT peer-reviewed.
    """
    data = await crossref_search_chemrxiv(
        query=params.query,
        rows=params.limit,
        sort=params.sort,
        date_from=params.date_from,
        date_to=params.date_to,
    )

    if not data or not data.get("message", {}).get("items"):
        return f"No ChemRxiv preprints found for '{params.query}'."

    items = data["message"]["items"]
    total = data["message"].get("total-results", len(items))

    # Enrich with OpenAlex abstracts
    for item in items:
        if not item.get("abstract") and item.get("DOI"):
            try:
                oa = await openalex_get_work(item["DOI"])
                if oa and oa.get("abstract_inverted_index"):
                    item["abstract"] = openalex_reconstruct_abstract(
                        oa["abstract_inverted_index"]
                    )
            except Exception:
                pass

    lines = [f"## ChemRxiv Preprints: '{params.query or 'recent'}'"]
    lines.append(f"Showing {len(items)} of {total:,} results ⚠️ *Not peer-reviewed*\n")
    for i, item in enumerate(items, 1):
        lines.append(_format_crossref_item_compact(item, index=i))
        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 10: get_chemrxiv_categories
# =============================================================================


@mcp.tool()
async def get_chemrxiv_categories() -> str:
    """List all ChemRxiv subject categories."""
    lines = ["## ChemRxiv Subject Categories\n"]
    for cid, name in sorted(CHEMRXIV_CATEGORIES.items(), key=lambda x: x[1]):
        lines.append(f"  - {name} (ID: {cid})")
    return "\n".join(lines)


# =============================================================================
# Tool 11: search_compound — PubChem + CAS compound lookup
# =============================================================================


class SearchCompoundInput(BaseModel):
    """Search for chemical compounds."""
    query: str = Field(
        description="Compound name, CAS number, SMILES, InChI, or molecular formula."
    )
    search_type: str = Field(
        default="auto",
        description=(
            "Search type: 'auto' (detect from query), 'name', 'cas', "
            "'smiles', 'formula'. Default: auto."
        ),
    )


@mcp.tool()
async def search_compound(params: SearchCompoundInput) -> str:
    """Search for chemical compounds by name, CAS number, SMILES, or formula.
    Returns structure information, identifiers, and basic properties.

    Data sources: PubChem (115M+ compounds) and CAS Common Chemistry (~500k).
    """
    query = params.query.strip()
    search_type = params.search_type.lower()

    # Auto-detect search type
    if search_type == "auto":
        import re
        if re.match(r"^\d{1,7}-\d{2}-\d$", query):
            search_type = "cas"
        elif any(c in query for c in "()=#@+\\[]") or query.startswith("["):
            search_type = "smiles"
        elif re.match(r"^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$", query) and len(query) < 30:
            search_type = "formula"
        else:
            search_type = "name"

    lines = [f"## Compound Search: {query}\n"]

    # --- CAS Common Chemistry ---
    cas_data = await cas_search(query)
    cas_results = []
    if cas_data and cas_data.get("results"):
        cas_results = cas_data["results"]
        lines.append("### CAS Common Chemistry")
        for r in cas_results[:5]:
            rn = r.get("rn", "")
            name = r.get("name", "")
            lines.append(f"  - **{name}** (CAS: {rn})")
        lines.append("")

        # Get details for first result
        if cas_results:
            first_cas = cas_results[0].get("rn", "")
            if first_cas:
                detail = await cas_detail(first_cas)
                if detail:
                    lines.append(f"### CAS Details: {detail.get('name', first_cas)}")
                    if detail.get("molecularFormula"):
                        lines.append(f"  **Formula:** {detail['molecularFormula']}")
                    if detail.get("molecularMass"):
                        lines.append(f"  **MW:** {detail['molecularMass']}")
                    if detail.get("inchi"):
                        lines.append(f"  **InChI:** {detail['inchi']}")
                    if detail.get("inchiKey"):
                        lines.append(f"  **InChIKey:** {detail['inchiKey']}")
                    if detail.get("smile"):
                        lines.append(f"  **SMILES:** {detail['smile']}")
                    if detail.get("canonicalSmile"):
                        lines.append(f"  **Canonical SMILES:** {detail['canonicalSmile']}")
                    # Experimental properties
                    props = detail.get("experimentalProperties", [])
                    if props:
                        lines.append("  **Experimental properties:**")
                        for p in props:
                            pname = p.get("name", "")
                            pval = p.get("property", "")
                            lines.append(f"    - {pname}: {pval}")
                    lines.append("")

    # --- PubChem ---
    pc_data = None
    if search_type == "smiles":
        pc_data = await pubchem_search_by_smiles(query)
    elif search_type == "formula":
        pc_data = await pubchem_search_by_formula(query)
    else:
        pc_data = await pubchem_search_by_name(query)

    if pc_data and pc_data.get("PC_Compounds"):
        compounds = pc_data["PC_Compounds"]
        lines.append("### PubChem")
        for comp in compounds[:3]:
            cid = comp.get("id", {}).get("id", {}).get("cid", "")
            if cid:
                # Get properties
                prop_data = await pubchem_get_properties(cid)
                if prop_data and prop_data.get("PropertyTable", {}).get("Properties"):
                    props = prop_data["PropertyTable"]["Properties"][0]
                    lines.append(f"  **PubChem CID:** {cid}")
                    if props.get("IUPACName"):
                        lines.append(f"  **IUPAC Name:** {props['IUPACName']}")
                    if props.get("MolecularFormula"):
                        lines.append(f"  **Formula:** {props['MolecularFormula']}")
                    if props.get("MolecularWeight"):
                        lines.append(f"  **MW:** {props['MolecularWeight']}")
                    if props.get("CanonicalSMILES"):
                        lines.append(f"  **SMILES:** {props['CanonicalSMILES']}")
                    if props.get("InChI"):
                        lines.append(f"  **InChI:** {props['InChI']}")
                    if props.get("InChIKey"):
                        lines.append(f"  **InChIKey:** {props['InChIKey']}")
                    if props.get("XLogP") is not None:
                        lines.append(f"  **XLogP:** {props['XLogP']}")
                    if props.get("TPSA") is not None:
                        lines.append(f"  **TPSA:** {props['TPSA']} Å²")
                    if props.get("HBondDonorCount") is not None:
                        lines.append(f"  **H-bond donors:** {props['HBondDonorCount']}")
                    if props.get("HBondAcceptorCount") is not None:
                        lines.append(f"  **H-bond acceptors:** {props['HBondAcceptorCount']}")
                    if props.get("RotatableBondCount") is not None:
                        lines.append(f"  **Rotatable bonds:** {props['RotatableBondCount']}")
                    lines.append(f"  **PubChem URL:** https://pubchem.ncbi.nlm.nih.gov/compound/{cid}")
                    lines.append("")

    if len(lines) <= 2:
        return f"No compounds found for '{query}'."

    return "\n".join(lines)


# =============================================================================
# Tool 12: get_compound_properties — Detailed PubChem properties
# =============================================================================


class GetCompoundPropsInput(BaseModel):
    """Get detailed molecular properties for a compound."""
    identifier: str = Field(
        description="PubChem CID (number) or compound name."
    )


@mcp.tool()
async def get_compound_properties(params: GetCompoundPropsInput) -> str:
    """Get detailed computed molecular properties from PubChem.
    Includes Lipinski rule-of-5 parameters, TPSA, XLogP, and all common names.
    """
    identifier = params.identifier.strip()

    # Resolve to CID if name given
    cid = None
    try:
        cid = int(identifier)
    except ValueError:
        # Search by name
        data = await pubchem_search_by_name(identifier)
        if data and data.get("PC_Compounds"):
            cid = data["PC_Compounds"][0].get("id", {}).get("id", {}).get("cid")

    if not cid:
        return f"Compound not found: '{identifier}'"

    # Fetch properties and synonyms in parallel
    props_task = pubchem_get_properties(cid)
    syn_task = pubchem_get_synonyms(cid)
    props_data, syn_data = await asyncio.gather(
        props_task, syn_task, return_exceptions=True
    )

    lines = [f"## Compound Properties (PubChem CID: {cid})\n"]

    if props_data and not isinstance(props_data, Exception):
        props_table = props_data.get("PropertyTable", {}).get("Properties", [])
        if props_table:
            p = props_table[0]
            if p.get("IUPACName"):
                lines.append(f"**IUPAC Name:** {p['IUPACName']}")
            if p.get("MolecularFormula"):
                lines.append(f"**Molecular Formula:** {p['MolecularFormula']}")
            if p.get("MolecularWeight"):
                lines.append(f"**Molecular Weight:** {p['MolecularWeight']} g/mol")
            if p.get("ExactMass"):
                lines.append(f"**Exact Mass:** {p['ExactMass']} Da")
            lines.append("")

            lines.append("### Structural Identifiers")
            if p.get("CanonicalSMILES"):
                lines.append(f"**Canonical SMILES:** `{p['CanonicalSMILES']}`")
            if p.get("IsomericSMILES"):
                lines.append(f"**Isomeric SMILES:** `{p['IsomericSMILES']}`")
            if p.get("InChI"):
                lines.append(f"**InChI:** `{p['InChI']}`")
            if p.get("InChIKey"):
                lines.append(f"**InChIKey:** `{p['InChIKey']}`")
            lines.append("")

            lines.append("### Drug-likeness (Lipinski Rule-of-5)")
            if p.get("MolecularWeight"):
                mw = float(p["MolecularWeight"])
                flag = "✅" if mw < 500 else "❌"
                lines.append(f"  {flag} MW: {mw:.1f} (< 500)")
            if p.get("XLogP") is not None:
                logp = p["XLogP"]
                flag = "✅" if logp < 5 else "❌"
                lines.append(f"  {flag} XLogP: {logp} (< 5)")
            if p.get("HBondDonorCount") is not None:
                hbd = p["HBondDonorCount"]
                flag = "✅" if hbd <= 5 else "❌"
                lines.append(f"  {flag} H-bond donors: {hbd} (≤ 5)")
            if p.get("HBondAcceptorCount") is not None:
                hba = p["HBondAcceptorCount"]
                flag = "✅" if hba <= 10 else "❌"
                lines.append(f"  {flag} H-bond acceptors: {hba} (≤ 10)")
            lines.append("")

            lines.append("### Additional Descriptors")
            if p.get("TPSA") is not None:
                lines.append(f"**TPSA:** {p['TPSA']} Å²")
            if p.get("Complexity") is not None:
                lines.append(f"**Complexity:** {p['Complexity']}")
            if p.get("RotatableBondCount") is not None:
                lines.append(f"**Rotatable Bonds:** {p['RotatableBondCount']}")
            if p.get("HeavyAtomCount") is not None:
                lines.append(f"**Heavy Atoms:** {p['HeavyAtomCount']}")
            lines.append("")

    # Synonyms / common names
    if syn_data and not isinstance(syn_data, Exception):
        syn_table = syn_data.get("InformationList", {}).get("Information", [])
        if syn_table and syn_table[0].get("Synonym"):
            names = syn_table[0]["Synonym"][:15]
            lines.append("### Common Names & Synonyms")
            lines.append(", ".join(names))
            if len(syn_table[0]["Synonym"]) > 15:
                lines.append(f"  ... (+{len(syn_table[0]['Synonym']) - 15} more)")
            lines.append("")

    lines.append(f"**PubChem:** https://pubchem.ncbi.nlm.nih.gov/compound/{cid}")

    return "\n".join(lines)


# =============================================================================
# Tool 13: search_web_of_science — Optional, credential-gated
# =============================================================================


class SearchWoSInput(BaseModel):
    """Search Web of Science (requires institutional API key)."""
    query: str = Field(
        description=(
            "Web of Science advanced search query. Examples: "
            "'TS=(asymmetric catalysis)', 'AU=(Maruoka)', "
            "'SO=(Nature Chemistry)', 'TS=(CRISPR) AND PY=(2023-2025)'. "
            "TS=topic, AU=author, SO=source, PY=year, DO=DOI."
        )
    )
    limit: int = Field(default=10, ge=1, le=50, description="Max results (1-50).")
    sort: str = Field(
        default="RS",
        description="Sort: RS=relevance, PY=year, TC=times cited, LD=load date."
    )


@mcp.tool()
async def search_web_of_science(params: SearchWoSInput) -> str:
    """Search Web of Science for peer-reviewed literature.

    REQUIRES: WOS_API_KEY environment variable (from developer.clarivate.com).
    Available free for institutions with Web of Science subscriptions.
    """
    if not wos_available():
        return (
            "⚠️ Web of Science API key not configured.\n\n"
            "To enable WoS search:\n"
            "1. Get an API key at https://developer.clarivate.com\n"
            "2. Set environment variable: WOS_API_KEY=your_key\n"
            "3. Restart the MCP server.\n\n"
            "Requires institutional Web of Science subscription."
        )

    data = await wos_search(
        query=params.query,
        limit=params.limit,
        sort_field=params.sort,
    )

    if not data:
        return f"Web of Science search failed for: {params.query}"

    # Parse WoS response format
    records = data.get("Data", {}).get("Records", {}).get("records", {}).get("REC", [])
    query_result = data.get("QueryResult", {})
    total = query_result.get("RecordsFound", 0)

    if not records:
        return f"No Web of Science results for: {params.query}"

    lines = [f"## Web of Science: {params.query}"]
    lines.append(f"Showing {len(records)} of {total:,} results\n")

    for i, rec in enumerate(records, 1):
        static = rec.get("static_data", {})
        summary = static.get("summary", {})
        titles = summary.get("titles", {}).get("title", [])
        title = next((t.get("content", "") for t in titles if t.get("type") == "item"), "Untitled")

        # Authors
        names_data = summary.get("names", {}).get("name", [])
        author_names = [n.get("display_name", "") for n in names_data[:5] if n.get("display_name")]

        # Source
        pub_info = summary.get("pub_info", {})
        source = next(
            (t.get("content", "") for t in titles if t.get("type") == "source"),
            ""
        )
        year = pub_info.get("pubyear", "")

        # Citation count
        dynamic = rec.get("dynamic_data", {})
        tc = dynamic.get("citation_related", {}).get("tc_list", {}).get("silo_tc", {})
        cite_count = tc.get("local_count", 0) if isinstance(tc, dict) else 0

        # DOI
        identifiers = static.get("item", {}).get("ids", {})
        doi = ""
        other_ids = static.get("dynamic_data", {}).get("cluster_related", {}).get("identifiers", {}).get("identifier", [])
        for oid in other_ids if isinstance(other_ids, list) else []:
            if oid.get("type") == "doi":
                doi = oid.get("value", "")

        lines.append(f"**{i}.** {title}")
        parts = []
        if year:
            parts.append(str(year))
        if source:
            parts.append(source)
        if cite_count:
            parts.append(f"{cite_count} citations")
        if parts:
            lines.append(f"  {' | '.join(parts)}")
        if author_names:
            lines.append(f"  Authors: {', '.join(author_names)}")
        if doi:
            lines.append(f"  DOI: {doi}")
        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 14: search_nist_webbook — Thermochemical data (scraping)
# =============================================================================


class SearchNISTInput(BaseModel):
    """Search NIST Chemistry WebBook for thermochemical data."""
    query: str = Field(
        description="Compound name, CAS number (e.g., '71-43-2'), or molecular formula."
    )
    search_type: str = Field(
        default="auto",
        description="Search type: 'auto' (detect), 'name', 'cas', 'formula'."
    )


@mcp.tool()
async def search_nist_webbook(params: SearchNISTInput) -> str:
    """Search the NIST Chemistry WebBook for thermochemical and physical property data.
    Returns melting/boiling points, enthalpies of formation/vaporization/fusion,
    heat capacities, entropies, and links to IR/MS/UV-Vis spectra.

    Data source: webbook.nist.gov (scraped — no official API).
    """
    import re as _re

    query = params.query.strip()
    search_type = params.search_type.lower()

    # Auto-detect
    if search_type == "auto":
        if _re.match(r"^\d{1,7}-\d{2}-\d$", query):
            search_type = "cas"
        elif _re.match(r"^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$", query) and len(query) < 25:
            search_type = "formula"
        else:
            search_type = "name"

    html = await nist_search(query, search_type)
    if not html:
        return f"Could not reach NIST WebBook for '{query}'."

    # Check if we hit a single compound page or search results
    if not nist_is_compound_page(html):
        results = nist_parse_search_results(html)
        if not results:
            return f"No NIST WebBook results for '{query}'."

        # Fetch the first compound
        first = results[0]
        from labmate_mcp.apis import nist_fetch
        html2 = await nist_fetch({"ID": first["nist_id"], "Mask": "FFF", "Units": "SI"})
        if html2 and nist_is_compound_page(html2):
            html = html2
        else:
            # Return search results list
            lines = [f"## NIST WebBook Search: '{query}'"]
            lines.append(f"Found {len(results)} compounds:\n")
            for i, r in enumerate(results[:20], 1):
                lines.append(f"  {i}. **{r['name']}** (NIST: {r['nist_id']})")
            return "\n".join(lines)

    info = nist_parse_compound(html)
    if not info:
        return f"Could not parse NIST data for '{query}'."

    lines = [f"## NIST WebBook: {info.get('name', query)}\n"]

    # Basic identifiers
    if info.get("formula"):
        lines.append(f"**Formula:** {info['formula']}")
    if info.get("molecular_weight"):
        lines.append(f"**Molecular Weight:** {info['molecular_weight']}")
    if info.get("cas_rn"):
        lines.append(f"**CAS:** {info['cas_rn']}")
    if info.get("inchi"):
        lines.append(f"**InChI:** `{info['inchi']}`")
    if info.get("inchi_key"):
        lines.append(f"**InChIKey:** `{info['inchi_key']}`")
    if info.get("other_names"):
        lines.append(f"**Other names:** {'; '.join(info['other_names'][:5])}")
    lines.append("")

    # Thermochemical data
    has_thermo = False
    thermo_fields = [
        ("delta_fH_gas_kJ_mol", "ΔfH° (gas)"),
        ("S_gas_J_mol_K", "S° (gas)"),
        ("Cp_gas_J_mol_K", "Cp (gas)"),
        ("boiling_point_K", "Boiling point"),
        ("melting_point_K", "Melting point"),
        ("delta_vapH_kJ_mol", "ΔvapH°"),
        ("delta_fusH_kJ_mol", "ΔfusH°"),
    ]
    for key, label in thermo_fields:
        if info.get(key):
            if not has_thermo:
                lines.append("### Thermochemical Data")
                has_thermo = True
            unit = ""
            if "kJ_mol" in key:
                unit = " kJ/mol"
            elif "J_mol_K" in key:
                unit = " J/(mol·K)"
            elif "_K" in key:
                unit = " K"
                # Also show Celsius
                try:
                    k_val = float(info[key])
                    c_val = k_val - 273.15
                    lines.append(f"  **{label}:** {info[key]} K ({c_val:.1f} °C)")
                    continue
                except ValueError:
                    pass
            lines.append(f"  **{label}:** {info[key]}{unit}")

    if has_thermo:
        lines.append("")

    # Available data
    avail = info.get("available_data", [])
    if avail:
        lines.append("### Available Data")
        nist_id = info.get("nist_id", "")
        for dtype in avail:
            url = f"https://webbook.nist.gov/cgi/cbook.cgi?ID={nist_id}&Units=SI"
            label = dtype.replace("_", " ").title()
            lines.append(f"  - {label}")
        if nist_id:
            lines.append(f"\n**Full NIST page:** https://webbook.nist.gov/cgi/cbook.cgi?ID={nist_id}&Units=SI&Mask=FFF")

    return "\n".join(lines)


# =============================================================================
# Tool 15: search_materials_project — Computational materials (credential-gated)
# =============================================================================


class SearchMaterialsInput(BaseModel):
    """Search Materials Project for inorganic materials."""
    formula: Optional[str] = Field(
        default=None,
        description="Chemical formula (e.g., 'LiFePO4', 'Fe2O3', 'SiO2'). Supports wildcards: 'Li*O*'."
    )
    elements: Optional[str] = Field(
        default=None,
        description="Comma-separated required elements (e.g., 'Li,Fe,O'). Returns all materials containing these."
    )
    band_gap_min: Optional[float] = Field(
        default=None, description="Minimum band gap in eV."
    )
    band_gap_max: Optional[float] = Field(
        default=None, description="Maximum band gap in eV."
    )
    limit: int = Field(default=10, ge=1, le=50, description="Max results.")


@mcp.tool()
async def search_materials_project(params: SearchMaterialsInput) -> str:
    """Search Materials Project for inorganic materials with DFT-calculated properties.
    Returns band gaps, formation energies, stability, and crystal symmetry.

    REQUIRES: MP_API_KEY environment variable (free at materialsproject.org).
    Contains 150,000+ computed materials.
    """
    if not mp_available():
        return (
            "⚠️ Materials Project API key not configured.\n\n"
            "To enable: set MP_API_KEY environment variable.\n"
            "Get a free key at https://materialsproject.org (register → dashboard → API key)."
        )

    elements = params.elements.split(",") if params.elements else None
    data = await mp_search(
        formula=params.formula,
        elements=[e.strip() for e in elements] if elements else None,
        band_gap_min=params.band_gap_min,
        band_gap_max=params.band_gap_max,
        limit=params.limit,
    )

    if not data or not data.get("data"):
        query_desc = params.formula or params.elements or "query"
        return f"No Materials Project results for '{query_desc}'."

    materials = data["data"]
    lines = [f"## Materials Project Results"]
    lines.append(f"Found {len(materials)} materials\n")

    for i, mat in enumerate(materials, 1):
        mp_id = mat.get("material_id", "")
        formula = mat.get("formula_pretty", "")
        bg = mat.get("band_gap")
        eform = mat.get("formation_energy_per_atom")
        ehull = mat.get("energy_above_hull")
        stable = mat.get("is_stable")
        sym = mat.get("symmetry", {})
        spg = sym.get("symbol", "") if isinstance(sym, dict) else ""
        crystal = sym.get("crystal_system", "") if isinstance(sym, dict) else ""

        lines.append(f"**{i}. {formula}** ({mp_id})")
        parts = []
        if bg is not None:
            parts.append(f"Band gap: {bg:.2f} eV")
        if eform is not None:
            parts.append(f"ΔHf: {eform:.3f} eV/atom")
        if ehull is not None:
            stability = "stable" if ehull == 0 else f"{ehull:.3f} eV above hull"
            parts.append(stability)
        if spg:
            parts.append(f"{spg}")
        if crystal:
            parts.append(crystal)
        if parts:
            lines.append(f"  {' | '.join(parts)}")
        lines.append(f"  https://materialsproject.org/materials/{mp_id}")
        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 16: predict_retrosynthesis — IBM RXN (credential-gated)
# =============================================================================


class RetrosynthesisInput(BaseModel):
    """AI-powered retrosynthetic analysis."""
    product_smiles: str = Field(
        description="SMILES of the target product molecule."
    )
    max_steps: int = Field(
        default=3, ge=1, le=10,
        description="Maximum retrosynthetic steps (1-10)."
    )


@mcp.tool()
async def predict_retrosynthesis(params: RetrosynthesisInput) -> str:
    """Predict retrosynthetic routes using IBM RXN for Chemistry.
    Uses transformer-based AI to plan multi-step synthesis pathways.

    REQUIRES: RXN_API_KEY environment variable (free at rxn.res.ibm.com).
    Rate limit: 5 requests/minute.
    """
    if not rxn_available():
        return (
            "⚠️ IBM RXN API key not configured.\n\n"
            "To enable: set RXN_API_KEY environment variable.\n"
            "Get a free key at https://rxn.res.ibm.com (register → API key)."
        )

    lines = [f"## Retrosynthetic Analysis"]
    lines.append(f"**Target:** `{params.product_smiles}`")
    lines.append(f"**Max steps:** {params.max_steps}")
    lines.append("\n*Running IBM RXN retrosynthesis (may take 30-120 seconds)...*\n")

    data = await rxn_retrosynthesis(params.product_smiles, params.max_steps)

    if not data:
        return lines[0] + "\n\n⚠️ Retrosynthesis timed out or failed. Try a simpler molecule or fewer steps."

    payload = data.get("payload", {})
    status = payload.get("status", data.get("status", "unknown"))

    if status not in ("SUCCESS", "success"):
        return lines[0] + f"\n\n⚠️ Retrosynthesis status: {status}"

    # Parse retrosynthetic trees
    retro_paths = payload.get("sequences", [])
    if not retro_paths:
        retro_paths = payload.get("retrosynthetic_paths", [])

    if not retro_paths:
        lines.append("No retrosynthetic routes found for this target.")
        return "\n".join(lines)

    lines.append(f"Found {len(retro_paths)} route(s):\n")

    for pi, path in enumerate(retro_paths[:5], 1):
        lines.append(f"### Route {pi}")
        reactions = path if isinstance(path, list) else path.get("reactions", [])
        for ri, rxn in enumerate(reactions, 1):
            if isinstance(rxn, dict):
                smiles = rxn.get("smiles", rxn.get("rxnSmiles", ""))
                confidence = rxn.get("confidence", "")
                conf_str = f" (confidence: {confidence:.2f})" if isinstance(confidence, (int, float)) else ""
                lines.append(f"  Step {ri}: `{smiles}`{conf_str}")
            else:
                lines.append(f"  Step {ri}: `{rxn}`")
        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 17: translate_compound_ids — UniChem cross-reference
# =============================================================================


class TranslateIDsInput(BaseModel):
    """Cross-reference compound identifiers across databases."""
    inchikey: str = Field(
        description="InChIKey of the compound (e.g., 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N' for aspirin)."
    )


@mcp.tool()
async def translate_compound_ids(params: TranslateIDsInput) -> str:
    """Translate a compound's InChIKey to identifiers in 40+ databases.
    Links PubChem CID ↔ ChEMBL ID ↔ DrugBank ID ↔ ZINC ID ↔ CAS RN and more.

    Uses EBI UniChem (no authentication required).
    """
    data = await unichem_lookup(params.inchikey.strip())

    if not data:
        return f"No UniChem results for InChIKey: {params.inchikey}"

    compounds = data.get("compounds", [])
    if not compounds:
        return f"InChIKey not found in UniChem: {params.inchikey}"

    lines = [f"## Cross-database Identifiers for {params.inchikey}\n"]

    # Source name mapping for common databases
    source_names = {
        "1": "ChEMBL", "2": "DrugBank", "3": "PDBe",
        "4": "IUPHAR", "7": "ChEBI", "10": "ZINC",
        "14": "KEGG", "17": "PharmGKB", "18": "HMDB",
        "20": "SelleckChem", "21": "PubChem (SID)",
        "22": "PubChem (CID)", "23": "MolPort",
        "24": "IBM Patent Data", "25": "NMR-Shift DB",
        "26": "ChemicalBook", "27": "DrugCentral",
        "28": "CAS Common Chemistry", "29": "Recon",
        "31": "BindingDB", "33": "EPA CompTox",
        "34": "BRENDA", "36": "Probes & Drugs",
        "37": "SwissLipids", "38": "Rhea",
        "39": "DP COSMOS", "41": "NIKKAJI",
    }

    for comp in compounds[:1]:  # First compound match
        sources = comp.get("sources", [])
        if sources:
            for src in sorted(sources, key=lambda s: s.get("sourceName", "")):
                src_id = str(src.get("sourceId", ""))
                src_name = source_names.get(src_id, src.get("sourceName", f"Source {src_id}"))
                compound_id = src.get("compoundId", "")
                url = src.get("url", "")
                url_str = f" — {url}" if url else ""
                lines.append(f"  - **{src_name}:** {compound_id}{url_str}")

    if len(lines) <= 2:
        lines.append("No cross-references found.")

    return "\n".join(lines)


# =============================================================================
# Tool 18: search_crystal_structures — COD
# =============================================================================


class SearchCrystalInput(BaseModel):
    """Search for crystal structures in the Crystallography Open Database."""
    formula: Optional[str] = Field(
        default=None,
        description="Chemical formula in Hill notation with spaces (e.g., 'C6 H6', 'Fe2 O3')."
    )
    elements: Optional[str] = Field(
        default=None,
        description="Comma-separated required elements (e.g., 'Fe,O')."
    )
    text: Optional[str] = Field(
        default=None,
        description="Free text search in compound names."
    )
    limit: int = Field(default=15, ge=1, le=50, description="Max results.")


@mcp.tool()
async def search_crystal_structures(params: SearchCrystalInput) -> str:
    """Search the Crystallography Open Database (COD) for experimental crystal structures.
    Returns space groups, cell parameters, and CIF file links.
    500,000+ structures under CC0 license. No authentication required.
    """
    elements = [e.strip() for e in params.elements.split(",")] if params.elements else None
    data = await cod_search(
        formula=params.formula,
        elements=elements,
        text=params.text,
        limit=params.limit,
    )

    if not data:
        query = params.formula or params.elements or params.text or ""
        return f"No COD crystal structures found for '{query}'."

    structures = data if isinstance(data, list) else []
    if not structures:
        return "No crystal structures found."

    lines = [f"## COD Crystal Structures"]
    lines.append(f"Found {len(structures)} structures\n")

    for i, s in enumerate(structures[:params.limit], 1):
        if isinstance(s, dict):
            cod_id = s.get("file", s.get("cod_id", ""))
            formula = s.get("formula", "")
            sg = s.get("sg", s.get("spacegroup", ""))
            a = s.get("a", "")
            b = s.get("b", "")
            c = s.get("c", "")
            title = s.get("title", "")

            lines.append(f"**{i}. COD {cod_id}** — {formula}")
            if title:
                lines.append(f"  {_truncate(title, 120)}")
            parts = []
            if sg:
                parts.append(f"Space group: {sg}")
            if a and b and c:
                parts.append(f"a={a}, b={b}, c={c} Å")
            if parts:
                lines.append(f"  {' | '.join(parts)}")
            lines.append(f"  CIF: https://www.crystallography.net/cod/{cod_id}.cif")
            lines.append("")
        else:
            # Simple list of IDs
            lines.append(f"  {i}. COD {s} — https://www.crystallography.net/cod/{s}.cif")

    return "\n".join(lines)


# =============================================================================
# Tool 19: search_toxicity — EPA CompTox (credential-gated)
# =============================================================================


class SearchToxInput(BaseModel):
    """Search EPA CompTox Dashboard for toxicity and environmental data."""
    query: str = Field(
        description="Chemical name, CAS number, or DTXSID identifier."
    )


@mcp.tool()
async def search_toxicity(params: SearchToxInput) -> str:
    """Search the EPA CompTox Dashboard for toxicity, environmental fate,
    and physicochemical property data. Contains 1.2M+ chemicals with
    ToxCast bioactivity screening and hazard classifications.

    REQUIRES: COMPTOX_API_KEY environment variable.
    Get a free key by emailing ccte_api@epa.gov.
    """
    if not comptox_available():
        return (
            "⚠️ EPA CompTox API key not configured.\n\n"
            "To enable: set COMPTOX_API_KEY environment variable.\n"
            "Request a free key by emailing ccte_api@epa.gov."
        )

    query = params.query.strip()
    data = await comptox_search(query)

    if not data:
        return f"No CompTox results for '{query}'."

    # Handle list vs single result
    chemicals = data if isinstance(data, list) else [data]
    if not chemicals:
        return f"No CompTox results for '{query}'."

    lines = [f"## EPA CompTox Dashboard: {query}\n"]

    for chem in chemicals[:3]:
        dtxsid = chem.get("dtxsid", "")
        name = chem.get("preferredName", chem.get("name", ""))
        cas = chem.get("casrn", "")
        mw = chem.get("molWeight", chem.get("averageMass", ""))
        formula = chem.get("molFormula", "")
        smiles = chem.get("smiles", "")

        lines.append(f"### {name}")
        if dtxsid:
            lines.append(f"**DTXSID:** {dtxsid}")
        if cas:
            lines.append(f"**CAS:** {cas}")
        if formula:
            lines.append(f"**Formula:** {formula}")
        if mw:
            lines.append(f"**MW:** {mw}")
        if smiles:
            lines.append(f"**SMILES:** `{smiles}`")
        lines.append("")

        # Fetch detailed properties if we have DTXSID
        if dtxsid:
            props = await comptox_get_properties(dtxsid)
            if props:
                prop_list = props if isinstance(props, list) else [props]
                if prop_list:
                    lines.append("**Physicochemical Properties:**")
                    for p in prop_list[:10]:
                        if isinstance(p, dict):
                            pname = p.get("propertyName", p.get("name", ""))
                            pval = p.get("value", p.get("propertyValue", ""))
                            punits = p.get("units", p.get("propertyUnits", ""))
                            if pname and pval:
                                lines.append(f"  - {pname}: {pval} {punits}")
                    lines.append("")

            # Hazard data
            hazard = await comptox_get_hazard(dtxsid)
            if hazard:
                haz_list = hazard if isinstance(hazard, list) else [hazard]
                if haz_list:
                    lines.append("**Hazard Data:**")
                    for h in haz_list[:8]:
                        if isinstance(h, dict):
                            htype = h.get("hazardType", h.get("type", ""))
                            hval = h.get("hazardValue", h.get("value", ""))
                            if htype and hval:
                                lines.append(f"  - {htype}: {hval}")
                    lines.append("")

        if dtxsid:
            lines.append(f"**Dashboard:** https://comptox.epa.gov/dashboard/chemical/details/{dtxsid}")
        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 20: search_mass_spectra — MassBank EU
# =============================================================================


class SearchMSInput(BaseModel):
    """Search MassBank for reference mass spectra."""
    compound_name: Optional[str] = Field(
        default=None, description="Compound name to search."
    )
    formula: Optional[str] = Field(
        default=None, description="Molecular formula."
    )
    inchikey: Optional[str] = Field(
        default=None, description="InChIKey identifier."
    )
    exact_mass_min: Optional[float] = Field(
        default=None, description="Minimum exact mass (Da)."
    )
    exact_mass_max: Optional[float] = Field(
        default=None, description="Maximum exact mass (Da)."
    )
    limit: int = Field(default=15, ge=1, le=50, description="Max results.")


@mcp.tool()
async def search_mass_spectra(params: SearchMSInput) -> str:
    """Search MassBank EU for reference mass spectra (MS/MS).
    120,000+ curated spectra from 18,500 compounds.
    Use for compound identification and spectral matching.
    No authentication required.
    """
    if not any([params.compound_name, params.formula, params.inchikey,
                params.exact_mass_min]):
        return "Please provide at least one search criterion (name, formula, InChIKey, or mass range)."

    data = await massbank_search(
        compound_name=params.compound_name,
        formula=params.formula,
        inchikey=params.inchikey,
        exact_mass_min=params.exact_mass_min,
        exact_mass_max=params.exact_mass_max,
        limit=params.limit,
    )

    if not data:
        query = params.compound_name or params.formula or params.inchikey or ""
        return f"No MassBank spectra found for '{query}'."

    records = data if isinstance(data, list) else []
    if not records:
        return "No mass spectra found."

    lines = ["## MassBank Mass Spectra"]
    lines.append(f"Found {len(records)} spectra\n")

    for i, rec in enumerate(records[:params.limit], 1):
        if isinstance(rec, dict):
            accession = rec.get("accession", "")
            name = rec.get("compound", {}).get("names", [""])[0] if isinstance(rec.get("compound"), dict) else rec.get("title", "")
            formula = rec.get("compound", {}).get("formula", "") if isinstance(rec.get("compound"), dict) else rec.get("formula", "")
            mass = rec.get("compound", {}).get("mass", "") if isinstance(rec.get("compound"), dict) else rec.get("exactMass", "")
            inst_type = rec.get("instrument_type", rec.get("ac", {}).get("instrument_type", "")) if isinstance(rec, dict) else ""
            ms_type = rec.get("ms_type", rec.get("ac", {}).get("mass_spectrometry", {}).get("ms_type", "")) if isinstance(rec, dict) else ""

            lines.append(f"**{i}. {accession}** — {name}")
            parts = []
            if formula:
                parts.append(formula)
            if mass:
                parts.append(f"m/z {mass}")
            if ms_type:
                parts.append(ms_type)
            if inst_type:
                parts.append(inst_type)
            if parts:
                lines.append(f"  {' | '.join(str(p) for p in parts)}")
            if accession:
                lines.append(f"  https://massbank.eu/MassBank/RecordDisplay?id={accession}")
            lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 21: search_binding_data — BindingDB
# =============================================================================


class SearchBindingInput(BaseModel):
    """Search BindingDB for protein-ligand binding affinities."""
    uniprot_id: Optional[str] = Field(
        default=None,
        description="UniProt ID of the protein target (e.g., 'P00918' for carbonic anhydrase II)."
    )
    smiles: Optional[str] = Field(
        default=None,
        description="SMILES of a compound to find targets for."
    )
    cutoff: float = Field(
        default=10000,
        description="Binding affinity cutoff in nM (for UniProt search) or similarity threshold 0-1 (for SMILES search)."
    )


@mcp.tool()
async def search_binding_data(params: SearchBindingInput) -> str:
    """Search BindingDB for quantitative binding affinities (Ki, IC50, Kd, EC50).
    3.2 million measurements including patent-derived data.
    No authentication required.

    Search by protein target (UniProt ID) or by compound (SMILES).
    """
    if not params.uniprot_id and not params.smiles:
        return "Please provide either a UniProt ID (target search) or SMILES (compound search)."

    if params.uniprot_id:
        data = await bindingdb_by_target(params.uniprot_id, int(params.cutoff))
        search_desc = f"UniProt:{params.uniprot_id} (cutoff: {params.cutoff} nM)"
    else:
        cutoff_sim = params.cutoff if params.cutoff <= 1 else 0.8
        data = await bindingdb_by_smiles(params.smiles, cutoff_sim)
        search_desc = f"SMILES:{params.smiles[:50]} (similarity ≥ {cutoff_sim})"

    if not data:
        return f"No BindingDB results for {search_desc}."

    entries = data if isinstance(data, list) else []
    if not entries:
        return f"No binding data found for {search_desc}."

    lines = [f"## BindingDB: {search_desc}"]
    lines.append(f"Found {len(entries)} binding measurements\n")

    for i, entry in enumerate(entries[:20], 1):
        if isinstance(entry, dict):
            # Extract key fields (names vary by endpoint)
            target = entry.get("target_name", entry.get("Target Name", ""))
            compound = entry.get("monomerid", entry.get("Ligand SMILES", entry.get("smiles", "")))
            ki = entry.get("ki_nm", entry.get("Ki (nM)", ""))
            ic50 = entry.get("ic50_nm", entry.get("IC50 (nM)", ""))
            kd = entry.get("kd_nm", entry.get("Kd (nM)", ""))
            ec50 = entry.get("ec50_nm", entry.get("EC50 (nM)", ""))

            lines.append(f"**{i}.**")
            if target:
                lines.append(f"  Target: {target}")
            if compound:
                lines.append(f"  Compound: `{_truncate(str(compound), 80)}`")

            affinities = []
            if ki:
                affinities.append(f"Ki = {ki} nM")
            if ic50:
                affinities.append(f"IC50 = {ic50} nM")
            if kd:
                affinities.append(f"Kd = {kd} nM")
            if ec50:
                affinities.append(f"EC50 = {ec50} nM")
            if affinities:
                lines.append(f"  {' | '.join(affinities)}")
            lines.append("")

    if len(entries) > 20:
        lines.append(f"*... and {len(entries) - 20} more results.*")

    return "\n".join(lines)


# =============================================================================
# Tool 22: generate_bibtex — DOI to BibTeX reference
# =============================================================================


class GenerateBibTeXInput(BaseModel):
    """Generate BibTeX entries from DOIs."""
    dois: str = Field(
        description=(
            "One or more DOIs, separated by commas or newlines. "
            "E.g., '10.1038/s41586-020-2649-2, 10.1021/jacs.3c01510'"
        )
    )


@mcp.tool()
async def generate_bibtex(params: GenerateBibTeXInput) -> str:
    """Generate BibTeX entries from DOIs for use in LaTeX/BibTeX.
    Uses Crossref content negotiation — works for any DOI in Crossref.
    Supports batch: pass multiple DOIs separated by commas.
    """
    import re as _re
    raw = params.dois.strip()
    # Parse DOIs from various formats
    dois = [
        d.strip().removeprefix("https://doi.org/").removeprefix("http://doi.org/")
        for d in _re.split(r"[,\n;]+", raw)
        if d.strip()
    ]

    if not dois:
        return "No valid DOIs provided."

    results = await crossref_get_bibtex_batch(dois)

    lines = [f"## BibTeX ({len(dois)} reference{'s' if len(dois) > 1 else ''})\n"]
    lines.append("```bibtex")

    success = 0
    for doi, bib in results:
        if bib:
            lines.append(bib)
            lines.append("")
            success += 1
        else:
            lines.append(f"% ERROR: Could not retrieve BibTeX for {doi}")
            lines.append("")

    lines.append("```")

    if success < len(dois):
        lines.append(f"\n⚠️ {len(dois) - success}/{len(dois)} DOIs failed — they may not be in Crossref.")

    return "\n".join(lines)


# =============================================================================
# Tool 23: profile_compound — Orchestrated compound report
# =============================================================================


class ProfileCompoundInput(BaseModel):
    """Get a comprehensive compound profile from all available sources."""
    query: str = Field(
        description="Compound name, CAS number, or SMILES string."
    )


@mcp.tool()
async def profile_compound(params: ProfileCompoundInput) -> str:
    """Generate a comprehensive compound profile by querying all available databases.
    Combines: PubChem (structure, properties) + CAS (registry) + NIST WebBook
    (thermochemistry) + UniChem (cross-references) + GHS (safety data).

    This is a 'one query, get everything' tool for complete compound characterization.
    """
    query = params.query.strip()
    lines = [f"## Compound Profile: {query}\n"]

    # Step 1: PubChem search to get identifiers
    from labmate_mcp.apis import (
        pubchem_search_by_name, pubchem_get_properties, pubchem_get_synonyms,
        cas_search as _cas_search, nist_search as _nist_search,
        nist_is_compound_page, nist_parse_compound as _nist_parse,
        unichem_lookup as _unichem, pubchem_get_ghs as _get_ghs,
        parse_ghs_data as _parse_ghs, gnps_classify_compound as _classify,
    )

    pc_data = await pubchem_search_by_name(query)
    cid = None
    smiles = None
    inchikey = None

    if pc_data:
        compounds = (
            pc_data.get("PC_Compounds", [])
            if "PC_Compounds" in pc_data
            else []
        )
        if compounds:
            comp = compounds[0]
            cid = comp.get("id", {}).get("id", {}).get("cid")

            # Extract SMILES and InChIKey from properties
            for prop in comp.get("props", []):
                urn = prop.get("urn", {})
                label = urn.get("label", "")
                val = prop.get("value", {})
                if label == "SMILES" and urn.get("name") == "Canonical":
                    smiles = val.get("sval")
                elif label == "InChIKey":
                    inchikey = val.get("sval")

    if cid:
        lines.append(f"**PubChem CID:** [{cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})")

    # Step 2: Properties from PubChem
    if cid:
        props = await pubchem_get_properties(cid)
        if props and "PropertyTable" in props:
            ptable = props["PropertyTable"].get("Properties", [{}])[0]
            mw = ptable.get("MolecularWeight")
            mf = ptable.get("MolecularFormula")
            iupac = ptable.get("IUPACName")
            canon_smiles = ptable.get("CanonicalSMILES")
            xlogp = ptable.get("XLogP")
            tpsa = ptable.get("TPSA")
            hbd = ptable.get("HBondDonorCount")
            hba = ptable.get("HBondAcceptorCount")
            rb = ptable.get("RotatableBondCount")
            ik = ptable.get("InChIKey")

            if ik:
                inchikey = ik
            if canon_smiles:
                smiles = canon_smiles

            lines.append("\n### Identity")
            if iupac:
                lines.append(f"**IUPAC Name:** {iupac}")
            if mf:
                lines.append(f"**Formula:** {mf}")
            if mw:
                lines.append(f"**Molecular Weight:** {mw}")
            if smiles:
                lines.append(f"**SMILES:** `{smiles}`")
            if inchikey:
                lines.append(f"**InChIKey:** `{inchikey}`")

            lines.append("\n### Molecular Properties")
            prop_parts = []
            if xlogp is not None:
                prop_parts.append(f"XLogP: {xlogp}")
            if tpsa is not None:
                prop_parts.append(f"TPSA: {tpsa} Å²")
            if hbd is not None:
                prop_parts.append(f"H-bond donors: {hbd}")
            if hba is not None:
                prop_parts.append(f"H-bond acceptors: {hba}")
            if rb is not None:
                prop_parts.append(f"Rotatable bonds: {rb}")

            # Lipinski Rule-of-5
            if mw and xlogp is not None and hbd is not None and hba is not None:
                violations = 0
                try:
                    if float(mw) > 500: violations += 1
                    if float(xlogp) > 5: violations += 1
                    if int(hbd) > 5: violations += 1
                    if int(hba) > 10: violations += 1
                except (ValueError, TypeError):
                    pass
                ro5 = "✅ Passes" if violations <= 1 else f"❌ {violations} violations"
                prop_parts.append(f"Lipinski Ro5: {ro5}")

            for pp in prop_parts:
                lines.append(f"  {pp}")

    # Step 3: CAS Registry
    cas_data = await _cas_search(query)
    if cas_data and "results" in cas_data:
        cas_results = cas_data["results"]
        if cas_results:
            cas_rn = cas_results[0].get("rn", "")
            cas_name = cas_results[0].get("name", "")
            if cas_rn:
                lines.append(f"\n### CAS Registry")
                lines.append(f"**CAS RN:** {cas_rn}")
                if cas_name:
                    lines.append(f"**CAS Name:** {cas_name}")

    # Step 4: NIST Thermochemistry
    nist_html = await _nist_search(query, "name")
    if nist_html and nist_is_compound_page(nist_html):
        nist_info = _nist_parse(nist_html)
        if nist_info:
            thermo_fields = [
                ("delta_fH_gas_kJ_mol", "ΔfH° (gas)", " kJ/mol"),
                ("S_gas_J_mol_K", "S° (gas)", " J/(mol·K)"),
                ("Cp_gas_J_mol_K", "Cp (gas)", " J/(mol·K)"),
                ("boiling_point_K", "Boiling point", " K"),
                ("melting_point_K", "Melting point", " K"),
                ("delta_vapH_kJ_mol", "ΔvapH°", " kJ/mol"),
                ("delta_fusH_kJ_mol", "ΔfusH°", " kJ/mol"),
            ]
            has_nist = False
            for key, label, unit in thermo_fields:
                if nist_info.get(key):
                    if not has_nist:
                        lines.append("\n### Thermochemical Data (NIST)")
                        has_nist = True
                    val_str = str(nist_info[key])
                    if "_K" in key:
                        try:
                            c = float(val_str) - 273.15
                            val_str = f"{val_str} K ({c:.1f} °C)"
                            unit = ""
                        except ValueError:
                            pass
                    lines.append(f"  {label}: {val_str}{unit}")

    # Step 5: GHS Safety Data
    if cid:
        ghs_raw = await _get_ghs(cid)
        if ghs_raw:
            ghs = _parse_ghs(ghs_raw)
            has_safety = (
                ghs.get("signal_word")
                or ghs.get("hazard_statements")
                or ghs.get("pictograms")
            )
            if has_safety:
                lines.append("\n### Safety (GHS)")
                if ghs["signal_word"]:
                    lines.append(f"  **Signal word:** {ghs['signal_word']}")
                if ghs["pictograms"]:
                    lines.append(f"  **Pictograms:** {', '.join(ghs['pictograms'][:6])}")
                if ghs["hazard_statements"]:
                    lines.append(f"  **H-statements:** {'; '.join(ghs['hazard_statements'][:8])}")
                if ghs["precautionary_statements"]:
                    lines.append(f"  **P-statements:** {'; '.join(ghs['precautionary_statements'][:6])}")

    # Step 6: Natural Product Classification (if we have SMILES)
    if smiles:
        np_data = await _classify(smiles)
        if np_data:
            pathway = np_data.get("pathway_results", [])
            superclass = np_data.get("superclass_results", [])
            np_class = np_data.get("class_results", [])
            isglycoside = np_data.get("isglycoside", False)

            if pathway or superclass or np_class:
                lines.append("\n### Natural Product Classification (NPClassifier)")
                if pathway:
                    lines.append(f"  Pathway: {', '.join(pathway)}")
                if superclass:
                    lines.append(f"  Superclass: {', '.join(superclass)}")
                if np_class:
                    lines.append(f"  Class: {', '.join(np_class)}")
                if isglycoside:
                    lines.append(f"  Glycoside: Yes")

    # Step 7: Cross-references (UniChem)
    if inchikey:
        uc_data = await _unichem(inchikey)
        if uc_data and uc_data.get("compounds"):
            compounds_uc = uc_data["compounds"]
            if compounds_uc:
                sources = compounds_uc[0].get("sources", [])
                if sources:
                    # Pick interesting ones
                    interesting = {
                        "1": "ChEMBL", "2": "DrugBank", "7": "ChEBI",
                        "10": "ZINC", "14": "KEGG", "18": "HMDB",
                        "22": "PubChem", "31": "BindingDB",
                        "33": "EPA CompTox",
                    }
                    found = []
                    for src in sources:
                        sid = str(src.get("sourceId", ""))
                        if sid in interesting:
                            found.append(
                                f"{interesting[sid]}: {src.get('compoundId', '')}"
                            )
                    if found:
                        lines.append("\n### Database Cross-References")
                        for f_entry in found:
                            lines.append(f"  {f_entry}")

    return "\n".join(lines)


# =============================================================================
# Tool 24: get_safety_data — GHS hazard classification
# =============================================================================


class GetSafetyInput(BaseModel):
    """Get GHS safety data for a compound."""
    compound_name: str = Field(
        description="Compound name (e.g., 'methanol', 'sulfuric acid')."
    )


@mcp.tool()
async def get_safety_data(params: GetSafetyInput) -> str:
    """Get GHS hazard classification for a chemical compound.
    Returns signal word, hazard pictograms, H-statements, and P-statements.
    Data from PubChem's GHS classification records.
    """
    from labmate_mcp.apis import (
        pubchem_search_by_name, pubchem_get_ghs as _get_ghs,
        parse_ghs_data as _parse_ghs,
    )

    # Resolve name to CID
    pc_data = await pubchem_search_by_name(params.compound_name.strip())
    if not pc_data:
        return f"Compound '{params.compound_name}' not found in PubChem."

    cid = None
    compounds = pc_data.get("PC_Compounds", [])
    if compounds:
        cid = compounds[0].get("id", {}).get("id", {}).get("cid")

    if not cid:
        return f"Could not resolve '{params.compound_name}' to a PubChem CID."

    ghs_raw = await _get_ghs(cid)
    if not ghs_raw:
        return f"No GHS data available for '{params.compound_name}' (CID: {cid})."

    ghs = _parse_ghs(ghs_raw)

    lines = [f"## GHS Safety Data: {params.compound_name}\n"]
    lines.append(f"**PubChem CID:** [{cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=GHS-Classification)")

    if ghs["signal_word"]:
        emoji = "🔴" if ghs["signal_word"].lower() == "danger" else "🟡"
        lines.append(f"\n**Signal Word:** {emoji} **{ghs['signal_word']}**")

    if ghs["pictograms"]:
        lines.append(f"\n**GHS Pictograms:**")
        for p in ghs["pictograms"]:
            lines.append(f"  - {p}")

    if ghs["hazard_statements"]:
        lines.append(f"\n**Hazard Statements (H):**")
        for h in ghs["hazard_statements"]:
            lines.append(f"  - {h}")

    if ghs["precautionary_statements"]:
        lines.append(f"\n**Precautionary Statements (P):**")
        for p in ghs["precautionary_statements"]:
            lines.append(f"  - {p}")

    if not any([ghs["signal_word"], ghs["pictograms"], ghs["hazard_statements"]]):
        lines.append("\nNo GHS classification data found. This may mean the compound is not classified as hazardous, or data is not yet available in PubChem.")

    return "\n".join(lines)


# =============================================================================
# Tool 25: search_protein_structures — RCSB PDB
# =============================================================================


class SearchPDBInput(BaseModel):
    """Search RCSB Protein Data Bank for structures."""
    query: str = Field(
        description="Search query (e.g., 'HIV protease', 'insulin receptor', 'kinase inhibitor')."
    )
    limit: int = Field(default=10, ge=1, le=50, description="Max results.")


@mcp.tool()
async def search_protein_structures(params: SearchPDBInput) -> str:
    """Search the RCSB Protein Data Bank for experimental protein/nucleic acid structures.
    230,000+ structures including X-ray, cryo-EM, and NMR.
    No authentication required.
    """
    data = await pdb_search(params.query, limit=params.limit)

    if not data:
        return f"No PDB structures found for '{params.query}'."

    total = data.get("total_count", 0)
    results = data.get("result_set", [])

    if not results:
        return f"No PDB structures found for '{params.query}'."

    lines = [f"## PDB Structures: '{params.query}'"]
    lines.append(f"Found {total} structures (showing {len(results)})\n")

    # Fetch details for each hit
    for i, result in enumerate(results[:params.limit], 1):
        pdb_id = result.get("identifier", "")
        score = result.get("score", 0)

        entry = await pdb_get_entry(pdb_id)
        if entry:
            struct = entry.get("struct", {})
            title = struct.get("title", "")
            cell = entry.get("cell", {})
            exptl = entry.get("exptl", [{}])
            method = exptl[0].get("method", "") if exptl else ""
            resolution = entry.get("rcsb_entry_info", {}).get("resolution_combined", [None])
            res_val = resolution[0] if resolution else None
            deposit_date = entry.get("rcsb_accession_info", {}).get("deposit_date", "")
            entity_count = entry.get("rcsb_entry_info", {}).get("polymer_entity_count", 0)

            lines.append(f"**{i}. [{pdb_id}](https://www.rcsb.org/structure/{pdb_id})**")
            if title:
                lines.append(f"  {_truncate(title, 120)}")

            parts = []
            if method:
                parts.append(method)
            if res_val:
                parts.append(f"{res_val:.2f} Å")
            if entity_count:
                parts.append(f"{entity_count} chain{'s' if entity_count > 1 else ''}")
            if deposit_date:
                parts.append(deposit_date[:10])
            if parts:
                lines.append(f"  {' | '.join(parts)}")
            lines.append("")
        else:
            lines.append(f"**{i}. [{pdb_id}](https://www.rcsb.org/structure/{pdb_id})**\n")

    return "\n".join(lines)


# =============================================================================
# Tool 26: get_protein_structure — PDB entry details + ligands
# =============================================================================


class GetPDBInput(BaseModel):
    """Get detailed information about a PDB entry."""
    pdb_id: str = Field(
        description="PDB ID (e.g., '1HHP', '7BZ5', '6LU7')."
    )


@mcp.tool()
async def get_protein_structure(params: GetPDBInput) -> str:
    """Get full details for a specific PDB structure including method, resolution,
    chains, ligands, and primary citation.
    """
    pdb_id = params.pdb_id.strip().upper()
    entry = await pdb_get_entry(pdb_id)

    if not entry:
        return f"PDB entry '{pdb_id}' not found."

    lines = [f"## PDB Structure: {pdb_id}\n"]
    lines.append(f"**URL:** https://www.rcsb.org/structure/{pdb_id}")

    # Title
    struct = entry.get("struct", {})
    title = struct.get("title", "")
    if title:
        lines.append(f"**Title:** {title}")

    # Method & Resolution
    exptl = entry.get("exptl", [{}])
    method = exptl[0].get("method", "") if exptl else ""
    res = entry.get("rcsb_entry_info", {}).get("resolution_combined", [None])
    res_val = res[0] if res else None

    if method:
        lines.append(f"**Method:** {method}")
    if res_val:
        lines.append(f"**Resolution:** {res_val:.2f} Å")

    # Cell parameters
    cell = entry.get("cell", {})
    if cell:
        a, b, c = cell.get("length_a"), cell.get("length_b"), cell.get("length_c")
        sg = entry.get("symmetry", {}).get("space_group_name_H_M", "")
        if a and b and c:
            lines.append(f"**Unit cell:** a={a}, b={b}, c={c} Å")
        if sg:
            lines.append(f"**Space group:** {sg}")

    # Dates
    accession = entry.get("rcsb_accession_info", {})
    deposit = accession.get("deposit_date", "")
    release = accession.get("initial_release_date", "")
    if deposit:
        lines.append(f"**Deposited:** {deposit[:10]}")
    if release:
        lines.append(f"**Released:** {release[:10]}")

    # Primary citation
    citation = entry.get("rcsb_primary_citation", {})
    if citation:
        ctitle = citation.get("title", "")
        journal = citation.get("journal_abbrev", "")
        year = citation.get("year")
        doi = citation.get("pdbx_database_id_DOI", "")
        authors = citation.get("rcsb_authors", [])
        if ctitle:
            lines.append(f"\n### Primary Citation")
            lines.append(f"  {ctitle}")
            parts = []
            if authors:
                auth_str = ", ".join(authors[:3])
                if len(authors) > 3:
                    auth_str += " et al."
                parts.append(auth_str)
            if journal:
                parts.append(journal)
            if year:
                parts.append(str(year))
            if parts:
                lines.append(f"  {' | '.join(parts)}")
            if doi:
                lines.append(f"  https://doi.org/{doi}")

    # Polymer entities (protein chains)
    info = entry.get("rcsb_entry_info", {})
    poly_count = info.get("polymer_entity_count", 0)
    nonpoly_count = info.get("nonpolymer_entity_count", 0)

    if poly_count:
        lines.append(f"\n### Chains ({poly_count} polymer entit{'ies' if poly_count > 1 else 'y'})")
        for eid in range(1, min(poly_count + 1, 6)):
            ent = await pdb_get_entity(pdb_id, eid)
            if ent:
                desc = ent.get("rcsb_polymer_entity", {}).get("pdbx_description", "")
                src = ent.get("rcsb_entity_source_organism", [{}])
                organism = src[0].get("scientific_name", "") if src else ""
                seq_len = ent.get("entity_poly", {}).get("rcsb_sample_sequence_length", "")
                parts = []
                if desc:
                    parts.append(desc)
                if organism:
                    parts.append(organism)
                if seq_len:
                    parts.append(f"{seq_len} residues")
                if parts:
                    lines.append(f"  Entity {eid}: {' | '.join(parts)}")

    # Ligands
    if nonpoly_count:
        ligands = await pdb_get_ligands(pdb_id)
        if ligands:
            lines.append(f"\n### Ligands ({len(ligands)} bound molecule{'s' if len(ligands) > 1 else ''})")
            for lig in ligands:
                nonpoly = lig.get("rcsb_nonpolymer_entity", {})
                comp_id = lig.get("pdbx_entity_nonpoly", {}).get("comp_id", "")
                name = nonpoly.get("pdbx_description", "")
                formula = nonpoly.get("formula_weight", "")
                parts = []
                if comp_id:
                    parts.append(f"**{comp_id}**")
                if name:
                    parts.append(name)
                if formula:
                    parts.append(f"MW: {formula}")
                if parts:
                    lines.append(f"  {' — '.join(parts)}")

    return "\n".join(lines)


# =============================================================================
# Tool 27: classify_natural_product — GNPS NPClassifier
# =============================================================================


class ClassifyNPInput(BaseModel):
    """Classify a compound as a natural product."""
    smiles: str = Field(
        description="SMILES string of the compound to classify."
    )


@mcp.tool()
async def classify_natural_product(params: ClassifyNPInput) -> str:
    """Classify a compound into natural product categories using GNPS NPClassifier.
    Returns biosynthetic pathway, superclass, and class predictions.
    Uses ML models trained on natural product structural features.
    No authentication required.
    """
    data = await gnps_classify_compound(params.smiles.strip())

    if not data:
        return f"NPClassifier could not process SMILES: `{params.smiles}`"

    lines = [f"## Natural Product Classification\n"]
    lines.append(f"**SMILES:** `{params.smiles}`\n")

    pathway = data.get("pathway_results", [])
    superclass = data.get("superclass_results", [])
    np_class = data.get("class_results", [])
    isglycoside = data.get("isglycoside", False)

    if not any([pathway, superclass, np_class]):
        lines.append("No natural product classification could be determined. This compound may be purely synthetic or outside the NPClassifier training domain.")
        return "\n".join(lines)

    if pathway:
        lines.append(f"**Biosynthetic Pathway:** {', '.join(pathway)}")
    if superclass:
        lines.append(f"**Superclass:** {', '.join(superclass)}")
    if np_class:
        lines.append(f"**Class:** {', '.join(np_class)}")
    lines.append(f"**Glycoside:** {'Yes' if isglycoside else 'No'}")

    # Add hierarchy visualization
    if pathway and superclass and np_class:
        lines.append(f"\n**Classification hierarchy:**")
        lines.append(f"  {' → '.join(pathway)} → {' → '.join(superclass)} → {' → '.join(np_class)}")

    return "\n".join(lines)


# =============================================================================
# Tool 28: get_journal_metrics — Journal impact via OpenAlex
# =============================================================================


class GetJournalInput(BaseModel):
    """Get journal metrics and publication statistics."""
    query: str = Field(
        description="Journal name (e.g., 'Nature Chemistry', 'JACS', 'Angewandte Chemie')."
    )


@mcp.tool()
async def get_journal_metrics(params: GetJournalInput) -> str:
    """Look up journal metrics including h-index, citation statistics,
    open access percentage, and publication counts.
    Data from OpenAlex (free, no authentication required).
    """
    data = await openalex_search_sources(params.query.strip())

    if not data or not data.get("results"):
        return f"No journal found matching '{params.query}'."

    lines = [f"## Journal Metrics: '{params.query}'\n"]

    for i, source in enumerate(data["results"][:3], 1):
        name = source.get("display_name", "")
        oa_id = source.get("id", "")
        issn_l = source.get("issn_l", "")
        issn_list = source.get("issn", [])
        publisher = source.get("host_organization_name", "")
        src_type = source.get("type", "")
        is_oa = source.get("is_oa", False)
        works_count = source.get("works_count", 0)
        cited_by = source.get("cited_by_count", 0)
        h_index = source.get("summary_stats", {}).get("h_index", "")
        i10_index = source.get("summary_stats", {}).get("i10_index", "")
        two_yr_mean = source.get("summary_stats", {}).get("2yr_mean_citedness", "")
        homepage = source.get("homepage_url", "")

        lines.append(f"### {i}. {name}")
        if publisher:
            lines.append(f"**Publisher:** {publisher}")
        if issn_l:
            lines.append(f"**ISSN:** {issn_l}")
        if src_type:
            lines.append(f"**Type:** {src_type}")
        if is_oa:
            lines.append(f"**Open Access:** Yes ✅")
        if homepage:
            lines.append(f"**Homepage:** {homepage}")

        lines.append("")

        # Metrics
        metrics = []
        if h_index:
            metrics.append(f"**h-index:** {h_index}")
        if i10_index:
            metrics.append(f"**i10-index:** {i10_index}")
        if two_yr_mean:
            metrics.append(f"**2-year mean citedness:** {two_yr_mean:.2f}")
        if works_count:
            metrics.append(f"**Total works:** {works_count:,}")
        if cited_by:
            metrics.append(f"**Total citations:** {cited_by:,}")

        if metrics:
            lines.append("**Metrics:**")
            for m in metrics:
                lines.append(f"  {m}")

        # Publication trend (last 5 years from counts_by_year)
        counts = source.get("counts_by_year", [])
        if counts:
            lines.append("\n**Recent publication volume:**")
            for cy in sorted(counts, key=lambda x: x.get("year", 0), reverse=True)[:5]:
                yr = cy.get("year", "")
                wc = cy.get("works_count", 0)
                cc = cy.get("cited_by_count", 0)
                lines.append(f"  {yr}: {wc:,} papers, {cc:,} citations")

        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Tool 29: predict_product — IBM RXN forward reaction prediction
# =============================================================================


class PredictProductInput(BaseModel):
    """Forward reaction prediction — predict products from reactants."""
    reactants_smiles: str = Field(
        description=(
            "Reactants as SMILES. Use '.' to separate multiple reactants, "
            "'>>' to separate reactants from empty product side. "
            "Example: 'CC(=O)O.OCC>>' (acetic acid + ethanol → ?)"
        )
    )


@mcp.tool()
async def predict_product(params: PredictProductInput) -> str:
    """Predict reaction product(s) from reactants using IBM RXN AI.
    The complement to retrosynthesis — give it starting materials and it predicts what you'll get.

    REQUIRES: RXN_API_KEY environment variable.
    """
    if not rxn_available():
        return (
            "⚠️ IBM RXN API key not configured.\n\n"
            "To enable: set RXN_API_KEY environment variable.\n"
            "Get a key at https://rxn.res.ibm.com (Individual/Team plan)."
        )

    smiles = params.reactants_smiles
    if ">>" not in smiles:
        smiles = smiles + ">>"

    lines = ["## Forward Reaction Prediction"]
    lines.append(f"**Reactants:** `{smiles}`")
    lines.append("\n*Running IBM RXN forward prediction...*\n")

    data = await rxn_predict_reaction(smiles)

    if not data:
        return lines[0] + "\n\n⚠️ Prediction timed out or failed."

    payload = data.get("payload", {})
    status = payload.get("status", data.get("status", "unknown"))

    if status not in ("SUCCESS", "success"):
        return lines[0] + f"\n\n⚠️ Prediction status: {status}"

    # Parse predictions
    attempts = payload.get("attempts", [])
    if not attempts:
        # Try alternative response structures
        result = payload.get("result", payload)
        if isinstance(result, dict):
            attempts = result.get("attempts", [])

    if not attempts:
        lines.append("No product predictions returned.")
        return "\n".join(lines)

    lines.append("### Predicted Product(s)\n")
    for i, attempt in enumerate(attempts[:5], 1):
        if isinstance(attempt, dict):
            smiles_out = attempt.get("smiles", attempt.get("productMolecule", {}).get("smiles", ""))
            confidence = attempt.get("confidence", "")
            conf_str = f" — confidence: {confidence:.2%}" if isinstance(confidence, (int, float)) else ""
            lines.append(f"{i}. `{smiles_out}`{conf_str}")
        elif isinstance(attempt, str):
            lines.append(f"{i}. `{attempt}`")

    return "\n".join(lines)


# =============================================================================
# Tool 30: text_to_procedure — IBM RXN paragraph to actions
# =============================================================================


class TextToProcedureInput(BaseModel):
    """Convert experimental text to structured procedure steps."""
    paragraph: str = Field(
        description=(
            "Free-text experimental procedure paragraph. Example: "
            "'To a solution of compound A (100 mg) in THF (5 mL) was added "
            "NaH (60 mg) at 0°C. The mixture was stirred for 2h, then quenched "
            "with water and extracted with EtOAc.'"
        )
    )


@mcp.tool()
async def text_to_procedure(params: TextToProcedureInput) -> str:
    """Extract machine-readable steps from experimental text using IBM RXN NLP.
    Converts free-text experimental paragraphs into structured action sequences
    (MAKESOLUTION, ADD, STIR, FILTER, CONCENTRATE, etc.).

    Useful for: digitizing procedures from papers, standardizing experimental sections,
    thesis writing, ensuring reproducibility.

    REQUIRES: RXN_API_KEY environment variable.
    """
    if not rxn_available():
        return (
            "⚠️ IBM RXN API key not configured.\n\n"
            "To enable: set RXN_API_KEY environment variable.\n"
            "Get a key at https://rxn.res.ibm.com (Individual/Team plan)."
        )

    lines = ["## Text → Procedure Extraction"]
    lines.append(f"**Input text:** {params.paragraph[:200]}{'...' if len(params.paragraph) > 200 else ''}")
    lines.append("")

    data = await rxn_paragraph_to_actions(params.paragraph)

    if not data:
        return lines[0] + "\n\n⚠️ Extraction failed or timed out."

    payload = data.get("payload", data)
    actions = payload.get("actions", payload.get("output", []))

    if not actions:
        lines.append("No structured actions could be extracted from this text.")
        return "\n".join(lines)

    lines.append(f"### Extracted Actions ({len(actions)} steps)\n")
    for i, action in enumerate(actions, 1):
        if isinstance(action, dict):
            action_name = action.get("name", action.get("action", "UNKNOWN"))
            content = action.get("content", "")
            duration = action.get("duration", "")
            temperature = action.get("temperature", "")
            details = []
            if content:
                details.append(str(content))
            if duration:
                details.append(f"duration: {duration}")
            if temperature:
                details.append(f"temp: {temperature}")
            detail_str = f" — {', '.join(details)}" if details else ""
            lines.append(f"{i}. **{action_name}**{detail_str}")
        elif isinstance(action, str):
            lines.append(f"{i}. {action}")

    return "\n".join(lines)


# =============================================================================
# Tool 31: predict_atom_mapping — IBM RXN atom-to-atom mapping
# =============================================================================


class AtomMappingInput(BaseModel):
    """Atom-to-atom mapping for a reaction."""
    rxn_smiles: str = Field(
        description=(
            "Full reaction SMILES with reactants>>products. "
            "Example: 'CC(=O)O.OCC>>CC(=O)OCC.O' (esterification)"
        )
    )


@mcp.tool()
async def predict_atom_mapping(params: AtomMappingInput) -> str:
    """Map atoms from reactants to products using IBM RXN AI (atom-mapping-2020 model).
    Shows which atom in the reactants becomes which atom in the products.

    Useful for: understanding reaction mechanisms, reaction database curation,
    verifying proposed mechanisms.

    REQUIRES: RXN_API_KEY environment variable.
    """
    if not rxn_available():
        return (
            "⚠️ IBM RXN API key not configured.\n\n"
            "To enable: set RXN_API_KEY environment variable.\n"
            "Get a key at https://rxn.res.ibm.com (Individual/Team plan)."
        )

    lines = ["## Atom-to-Atom Mapping"]
    lines.append(f"**Reaction:** `{params.rxn_smiles}`")
    lines.append("\n*Running atom mapping model...*\n")

    data = await rxn_predict_atom_mapping(params.rxn_smiles)

    if not data:
        return lines[0] + "\n\n⚠️ Atom mapping timed out or failed."

    payload = data.get("payload", {})
    status = payload.get("status", data.get("status", "unknown"))

    if status not in ("SUCCESS", "success"):
        return lines[0] + f"\n\n⚠️ Status: {status}"

    # Parse atom mapping result
    attempts = payload.get("attempts", [])
    if not attempts:
        result = payload.get("result", payload)
        if isinstance(result, dict):
            attempts = result.get("attempts", [])

    if attempts:
        lines.append("### Mapped Reaction(s)\n")
        for i, attempt in enumerate(attempts[:3], 1):
            if isinstance(attempt, dict):
                mapped = attempt.get("smiles", attempt.get("mappedSmiles", ""))
                confidence = attempt.get("confidence", "")
                conf_str = f" (confidence: {confidence:.2%})" if isinstance(confidence, (int, float)) else ""
                lines.append(f"{i}. `{mapped}`{conf_str}")
            elif isinstance(attempt, str):
                lines.append(f"{i}. `{attempt}`")
    else:
        # Try to get mapping from top-level payload
        mapped_rxn = payload.get("mappedReactionSmiles", payload.get("mapped_smiles", ""))
        if mapped_rxn:
            lines.append(f"### Mapped Reaction\n\n`{mapped_rxn}`")
        else:
            lines.append("No atom mapping returned. Check that your reaction SMILES is valid.")

    lines.append("\n*Atom map numbers in the SMILES show corresponding atoms across reactants → products.*")

    return "\n".join(lines)


# =============================================================================
# Tool 32: plan_synthesis — IBM RXN synthesis planning (from retrosynthesis)
# =============================================================================


class SynthesisPlanInput(BaseModel):
    """Generate step-by-step synthesis procedure from a retrosynthesis result."""
    prediction_id: str = Field(
        description=(
            "Retrosynthesis prediction ID from a previous predict_retrosynthesis call. "
            "Found in the retrosynthesis result payload."
        )
    )
    sequence_index: int = Field(
        default=0, ge=0,
        description="Which retrosynthetic route to use (0 = first/best route)."
    )


@mcp.tool()
async def plan_synthesis(params: SynthesisPlanInput) -> str:
    """Turn a retrosynthetic route into a step-by-step synthesis procedure with actions.
    Chain this after predict_retrosynthesis to get a full synthesis plan.

    Workflow: predict_retrosynthesis → plan_synthesis → detailed procedure

    REQUIRES: RXN_API_KEY environment variable.
    """
    if not rxn_available():
        return (
            "⚠️ IBM RXN API key not configured.\n\n"
            "To enable: set RXN_API_KEY environment variable.\n"
            "Get a key at https://rxn.res.ibm.com (Individual/Team plan)."
        )

    lines = ["## Synthesis Plan"]
    lines.append(f"**Prediction ID:** `{params.prediction_id}`")
    lines.append(f"**Route:** #{params.sequence_index + 1}")
    lines.append("\n*Generating synthesis procedure (may take 30-90 seconds)...*\n")

    data = await rxn_synthesis_plan(params.prediction_id, params.sequence_index)

    if not data:
        return lines[0] + "\n\n⚠️ Synthesis planning failed."

    if data.get("error"):
        return lines[0] + f"\n\n⚠️ {data['error']}"

    lines.append(f"Using route {params.sequence_index + 1} of {data.get('total_sequences', '?')}")
    lines.append(f"**Synthesis ID:** `{data.get('synthesis_id', 'N/A')}`\n")

    # Parse procedure
    procedure = data.get("procedure")
    if procedure and isinstance(procedure, dict):
        payload = procedure.get("payload", procedure)
        steps = payload.get("steps", payload.get("actions", []))
        if steps:
            lines.append(f"### Procedure ({len(steps)} steps)\n")
            for i, step in enumerate(steps, 1):
                if isinstance(step, dict):
                    action = step.get("name", step.get("action", "Step"))
                    content = step.get("content", step.get("description", ""))
                    product = step.get("product", "")
                    lines.append(f"**Step {i}: {action}**")
                    if content:
                        lines.append(f"  {content}")
                    if product:
                        lines.append(f"  → Product: `{product}`")
                    lines.append("")
                elif isinstance(step, str):
                    lines.append(f"**Step {i}:** {step}\n")
        else:
            lines.append("Procedure generated but no detailed steps available.")
    else:
        # Fall back to plan data
        plan = data.get("plan", {})
        if plan:
            payload = plan.get("payload", plan)
            lines.append("### Synthesis Plan\n")
            tree = payload.get("sequences", payload.get("tree", []))
            if isinstance(tree, list):
                for item in tree:
                    if isinstance(item, dict):
                        rxn = item.get("smiles", item.get("rxnSmiles", ""))
                        if rxn:
                            lines.append(f"- `{rxn}`")
            elif isinstance(tree, dict):
                lines.append(f"```\n{tree}\n```")
        else:
            lines.append("Synthesis plan created but procedure details not available.")

    return "\n".join(lines)


# =============================================================================
# Tool 33: predict_pka — Rowan Science pKa prediction
# =============================================================================

_ROWAN_NOT_CONFIGURED = (
    "⚠️ Rowan Science not configured.\n\n"
    "To enable computational chemistry tools:\n"
    "1. `pip install rowan-python`\n"
    "2. Set ROWAN_API_KEY environment variable\n"
    "3. Get an API key at https://labs.rowansci.com/account/api-keys\n\n"
    "Note: Rowan workflows use credits. Check your balance at labs.rowansci.com."
)


def _rowan_format_header(title: str, smiles: str, credits: float | None = None) -> list[str]:
    """Common header for Rowan tool output."""
    lines = [f"## {title}"]
    lines.append(f"**Molecule:** `{smiles}`")
    lines.append(f"**Engine:** Rowan Science (cloud quantum chemistry)")
    if credits is not None:
        lines.append(f"**Credits used:** {credits:.2f}")
    return lines


class PredictPkaInput(BaseModel):
    """Predict acid/base dissociation constants."""
    smiles: str = Field(
        description="SMILES of the molecule (e.g., 'c1ccccc1O' for phenol)."
    )
    pka_range: Optional[str] = Field(
        default="2,12",
        description="pKa range to search as 'min,max' (default: '2,12')."
    )


@mcp.tool()
async def predict_pka(params: PredictPkaInput) -> str:
    """Predict pKa values using Rowan Science quantum chemistry.
    Computes acid and base dissociation constants for any molecule from SMILES.

    Uses AIMNet2 neural network potential with careful mode for accuracy.

    REQUIRES: rowan-python package + ROWAN_API_KEY (uses credits).
    """
    if not rowan_available():
        return _ROWAN_NOT_CONFIGURED

    # Parse pKa range
    try:
        parts = params.pka_range.split(",")
        pka_min, pka_max = int(parts[0].strip()), int(parts[1].strip())
    except Exception:
        pka_min, pka_max = 2, 12

    lines = _rowan_format_header("pKa Prediction", params.smiles)
    lines.append(f"**Range:** {pka_min}–{pka_max}")
    lines.append("\n*Computing pKa values (typically 10-60 seconds)...*\n")

    data = await rowan_predict_pka(params.smiles, pka_range=(pka_min, pka_max))

    if not data:
        return lines[0] + "\n\n⚠️ pKa computation failed."
    if data.get("error"):
        return lines[0] + f"\n\n⚠️ {data['error']}"

    lines[2] = f"**Engine:** Rowan Science — credits used: {data.get('credits_charged', '?')}"

    result = data.get("data", {})
    if not result:
        lines.append("Computation completed but no pKa data returned.")
        return "\n".join(lines)

    # Format pKa results
    lines.append("### Results\n")

    strongest_acid = result.get("strongest_acid")
    strongest_base = result.get("strongest_base")
    if strongest_acid is not None:
        lines.append(f"**Strongest acid pKa:** {strongest_acid:.2f}")
    if strongest_base is not None:
        lines.append(f"**Strongest base pKa:** {strongest_base:.2f}")

    # Individual sites
    sites = result.get("pka_sites", result.get("sites", result.get("results", [])))
    if isinstance(sites, list) and sites:
        lines.append(f"\n### Individual Sites ({len(sites)})\n")
        for site in sites:
            if isinstance(site, dict):
                pka = site.get("pka", site.get("value", ""))
                atom = site.get("atom_index", site.get("atom", ""))
                stype = site.get("type", site.get("site_type", ""))
                pka_str = f"{pka:.2f}" if isinstance(pka, (int, float)) else str(pka)
                lines.append(f"- pKa = {pka_str} (atom {atom}, {stype})")

    # Also try to get any other useful data
    for key in ["microstate_pkas", "protonation_states", "charge_states"]:
        if key in result and result[key]:
            lines.append(f"\n**{key.replace('_', ' ').title()}:** {result[key]}")

    return "\n".join(lines)


# =============================================================================
# Tool 34: predict_solubility — Rowan Science solubility prediction
# =============================================================================


class PredictSolubilityInput(BaseModel):
    """Predict aqueous or organic solvent solubility."""
    smiles: str = Field(
        description="SMILES of the molecule."
    )
    method: Optional[str] = Field(
        default="fastsolv",
        description="Prediction method: 'fastsolv' (default, fast ML), 'kingfisher', or 'esol'."
    )
    solvents: Optional[str] = Field(
        default=None,
        description="Comma-separated solvents (e.g., 'water,ethanol'). Defaults to water."
    )


@mcp.tool()
async def predict_solubility(params: PredictSolubilityInput) -> str:
    """Predict molecular solubility using Rowan Science.
    Supports aqueous and organic solvents using ML models.

    REQUIRES: rowan-python package + ROWAN_API_KEY (uses credits).
    """
    if not rowan_available():
        return _ROWAN_NOT_CONFIGURED

    solvents = [s.strip() for s in params.solvents.split(",")] if params.solvents else None

    lines = _rowan_format_header("Solubility Prediction", params.smiles)
    lines.append(f"**Method:** {params.method}")
    if solvents:
        lines.append(f"**Solvents:** {', '.join(solvents)}")
    lines.append("\n*Computing solubility...*\n")

    data = await rowan_predict_solubility(
        params.smiles,
        method=params.method,
        solvents=solvents,
    )

    if not data:
        return lines[0] + "\n\n⚠️ Solubility computation failed."
    if data.get("error"):
        return lines[0] + f"\n\n⚠️ {data['error']}"

    result = data.get("data", {})
    if not result:
        lines.append("Computation completed but no solubility data returned.")
        return "\n".join(lines)

    lines.append("### Results\n")
    if isinstance(result, dict):
        for key, value in result.items():
            if key.startswith("_"):
                continue
            if isinstance(value, (int, float)):
                lines.append(f"**{key.replace('_', ' ').title()}:** {value:.4f}")
            elif isinstance(value, list):
                lines.append(f"**{key.replace('_', ' ').title()}:**")
                for item in value:
                    lines.append(f"  - {item}")
            elif value is not None:
                lines.append(f"**{key.replace('_', ' ').title()}:** {value}")

    return "\n".join(lines)


# =============================================================================
# Tool 35: predict_admet — Rowan Science ADMET prediction
# =============================================================================


class PredictAdmetInput(BaseModel):
    """Predict ADMET (absorption, distribution, metabolism, excretion, toxicity) properties."""
    smiles: str = Field(
        description="SMILES of the molecule."
    )


@mcp.tool()
async def predict_admet(params: PredictAdmetInput) -> str:
    """Predict ADMET properties using Rowan Science.
    Computes absorption, distribution, metabolism, excretion, and toxicity predictions.

    Returns drug-likeness assessments useful for medicinal chemistry and pharmacology.

    REQUIRES: rowan-python package + ROWAN_API_KEY (uses credits).
    """
    if not rowan_available():
        return _ROWAN_NOT_CONFIGURED

    lines = _rowan_format_header("ADMET Prediction", params.smiles)
    lines.append("\n*Computing ADMET properties...*\n")

    data = await rowan_predict_admet(params.smiles)

    if not data:
        return lines[0] + "\n\n⚠️ ADMET computation failed."
    if data.get("error"):
        return lines[0] + f"\n\n⚠️ {data['error']}"

    result = data.get("data", {})
    if not result:
        lines.append("Computation completed but no ADMET data returned.")
        return "\n".join(lines)

    lines.append("### ADMET Properties\n")

    # Format ADMET categories
    categories = {
        "absorption": [], "distribution": [], "metabolism": [],
        "excretion": [], "toxicity": [], "other": [],
    }

    for key, value in result.items():
        if key.startswith("_"):
            continue
        key_lower = key.lower()
        placed = False
        for cat in ["absorption", "distribution", "metabolism", "excretion", "toxicity"]:
            if cat[:4] in key_lower:
                categories[cat].append((key, value))
                placed = True
                break
        if not placed:
            categories["other"].append((key, value))

    # If no categorization worked, just list everything
    has_categories = any(v for k, v in categories.items() if k != "other")
    if has_categories:
        for cat_name, props in categories.items():
            if props:
                lines.append(f"**{cat_name.title()}**")
                for k, v in props:
                    label = k.replace("_", " ").title()
                    if isinstance(v, float):
                        lines.append(f"  {label}: {v:.4f}")
                    else:
                        lines.append(f"  {label}: {v}")
                lines.append("")
    else:
        for key, value in result.items():
            if key.startswith("_"):
                continue
            label = key.replace("_", " ").title()
            if isinstance(value, float):
                lines.append(f"- **{label}:** {value:.4f}")
            elif isinstance(value, dict):
                lines.append(f"- **{label}:**")
                for k2, v2 in value.items():
                    lines.append(f"    {k2}: {v2}")
            else:
                lines.append(f"- **{label}:** {value}")

    return "\n".join(lines)


# =============================================================================
# Tool 36: search_tautomers — Rowan Science tautomer enumeration
# =============================================================================


class SearchTautomersInput(BaseModel):
    """Enumerate and rank tautomers of a molecule."""
    smiles: str = Field(
        description="SMILES of the molecule (e.g., 'C1=CC(=O)NC=C1' for 2-pyridinone)."
    )


@mcp.tool()
async def search_tautomers(params: SearchTautomersInput) -> str:
    """Enumerate and rank tautomers using Rowan Science.
    Uses quantum chemistry to identify all relevant tautomeric forms and rank them by energy.

    REQUIRES: rowan-python package + ROWAN_API_KEY (uses credits).
    """
    if not rowan_available():
        return _ROWAN_NOT_CONFIGURED

    lines = _rowan_format_header("Tautomer Search", params.smiles)
    lines.append("\n*Enumerating and ranking tautomers...*\n")

    data = await rowan_search_tautomers(params.smiles)

    if not data:
        return lines[0] + "\n\n⚠️ Tautomer search failed."
    if data.get("error"):
        return lines[0] + f"\n\n⚠️ {data['error']}"

    result = data.get("data", {})
    if not result:
        lines.append("Computation completed but no tautomer data returned.")
        return "\n".join(lines)

    lines.append("### Tautomers\n")

    tautomers = result.get("tautomers", result.get("conformers", []))
    if isinstance(tautomers, list) and tautomers:
        for i, taut in enumerate(tautomers, 1):
            if isinstance(taut, dict):
                smiles = taut.get("smiles", taut.get("canonical_smiles", ""))
                energy = taut.get("energy", taut.get("relative_energy", ""))
                pop = taut.get("population", taut.get("boltzmann_weight", ""))
                lines.append(f"**Tautomer {i}:** `{smiles}`")
                if isinstance(energy, (int, float)):
                    lines.append(f"  Relative energy: {energy:.2f} kcal/mol")
                if isinstance(pop, (int, float)):
                    lines.append(f"  Population: {pop:.1%}")
                lines.append("")
            elif isinstance(taut, str):
                lines.append(f"{i}. `{taut}`")
    else:
        # Generic data dump
        for key, value in result.items():
            if key.startswith("_"):
                continue
            label = key.replace("_", " ").title()
            if isinstance(value, (int, float)):
                lines.append(f"**{label}:** {value:.4f}")
            elif isinstance(value, list) and len(value) <= 10:
                lines.append(f"**{label}:**")
                for item in value:
                    lines.append(f"  - {item}")
            elif value is not None:
                lines.append(f"**{label}:** {value}")

    return "\n".join(lines)


# =============================================================================
# Tool 37: compute_descriptors — Rowan Science molecular descriptors
# =============================================================================


class ComputeDescriptorsInput(BaseModel):
    """Compute molecular descriptors from structure."""
    smiles: str = Field(
        description="SMILES of the molecule."
    )


@mcp.tool()
async def compute_descriptors(params: ComputeDescriptorsInput) -> str:
    """Compute molecular descriptors using Rowan Science.
    Returns cheminformatics descriptors computed from the molecular structure.

    REQUIRES: rowan-python package + ROWAN_API_KEY (uses credits).
    """
    if not rowan_available():
        return _ROWAN_NOT_CONFIGURED

    lines = _rowan_format_header("Molecular Descriptors", params.smiles)
    lines.append("\n*Computing descriptors...*\n")

    data = await rowan_compute_descriptors(params.smiles)

    if not data:
        return lines[0] + "\n\n⚠️ Descriptor computation failed."
    if data.get("error"):
        return lines[0] + f"\n\n⚠️ {data['error']}"

    result = data.get("data", {})
    if not result:
        lines.append("Computation completed but no descriptor data returned.")
        return "\n".join(lines)

    lines.append("### Computed Descriptors\n")

    # Try to group descriptors by category
    for key, value in sorted(result.items()):
        if key.startswith("_"):
            continue
        label = key.replace("_", " ").title()
        if isinstance(value, float):
            lines.append(f"- **{label}:** {value:.4f}")
        elif isinstance(value, int):
            lines.append(f"- **{label}:** {value}")
        elif isinstance(value, dict):
            lines.append(f"- **{label}:**")
            for k2, v2 in value.items():
                if isinstance(v2, float):
                    lines.append(f"    {k2}: {v2:.4f}")
                else:
                    lines.append(f"    {k2}: {v2}")
        elif value is not None:
            lines.append(f"- **{label}:** {value}")

    return "\n".join(lines)


# =============================================================================
# Tool 38: predict_nmr — Rowan Science NMR prediction
# =============================================================================


class PredictNmrInput(BaseModel):
    """Predict NMR chemical shifts."""
    smiles: str = Field(
        description="SMILES of the molecule."
    )
    solvent: Optional[str] = Field(
        default="chloroform",
        description="NMR solvent (default: 'chloroform'). Common: 'chloroform', 'dmso', 'water', 'methanol'."
    )


@mcp.tool()
async def predict_nmr(params: PredictNmrInput) -> str:
    """Predict ¹H and ¹³C NMR chemical shifts using Rowan Science.
    Uses conformer search + neural network potentials to predict NMR spectra.

    Useful for: structure verification, spectral assignment, thesis writing.

    REQUIRES: rowan-python package + ROWAN_API_KEY (uses credits).
    """
    if not rowan_available():
        return _ROWAN_NOT_CONFIGURED

    lines = _rowan_format_header("NMR Prediction", params.smiles)
    lines.append(f"**Solvent:** {params.solvent}")
    lines.append("\n*Computing NMR shifts (conformer search + prediction, may take 1-3 minutes)...*\n")

    data = await rowan_predict_nmr(params.smiles, solvent=params.solvent)

    if not data:
        return lines[0] + "\n\n⚠️ NMR prediction failed."
    if data.get("error"):
        return lines[0] + f"\n\n⚠️ {data['error']}"

    result = data.get("data", {})
    if not result:
        lines.append("Computation completed but no NMR data returned.")
        return "\n".join(lines)

    lines.append("### Predicted NMR Shifts\n")

    # Try common result structures
    h_shifts = result.get("h_shifts", result.get("proton_shifts", result.get("1H", [])))
    c_shifts = result.get("c_shifts", result.get("carbon_shifts", result.get("13C", [])))

    if h_shifts:
        lines.append("**¹H Chemical Shifts (ppm)**\n")
        if isinstance(h_shifts, list):
            for shift in h_shifts:
                if isinstance(shift, dict):
                    atom = shift.get("atom_index", shift.get("atom", ""))
                    val = shift.get("shift", shift.get("value", ""))
                    val_str = f"{val:.2f}" if isinstance(val, (int, float)) else str(val)
                    lines.append(f"  H{atom}: {val_str} ppm")
                elif isinstance(shift, (int, float)):
                    lines.append(f"  {shift:.2f} ppm")
        lines.append("")

    if c_shifts:
        lines.append("**¹³C Chemical Shifts (ppm)**\n")
        if isinstance(c_shifts, list):
            for shift in c_shifts:
                if isinstance(shift, dict):
                    atom = shift.get("atom_index", shift.get("atom", ""))
                    val = shift.get("shift", shift.get("value", ""))
                    val_str = f"{val:.2f}" if isinstance(val, (int, float)) else str(val)
                    lines.append(f"  C{atom}: {val_str} ppm")
                elif isinstance(shift, (int, float)):
                    lines.append(f"  {shift:.2f} ppm")
        lines.append("")

    # If the result structure is different, just dump it nicely
    if not h_shifts and not c_shifts:
        for key, value in result.items():
            if key.startswith("_"):
                continue
            label = key.replace("_", " ").title()
            if isinstance(value, list) and value:
                lines.append(f"**{label}:**")
                for item in value[:50]:
                    if isinstance(item, dict):
                        lines.append(f"  {item}")
                    elif isinstance(item, (int, float)):
                        lines.append(f"  {item:.2f}")
                    else:
                        lines.append(f"  {item}")
            elif isinstance(value, (int, float)):
                lines.append(f"**{label}:** {value:.4f}")
            elif value is not None:
                lines.append(f"**{label}:** {value}")

    return "\n".join(lines)


# =============================================================================
# Bench Chemistry — Calculators (5 tools, pure computation)
# =============================================================================


class CalcMolarityInput(BaseModel):
    """Calculate molarity, moles, mass, or volume — provide any 2–3 knowns, solve for unknowns."""
    concentration_M: float | None = Field(None, description="Concentration in mol/L (M)")
    mass_grams: float | None = Field(None, description="Mass in grams")
    moles: float | None = Field(None, description="Amount in moles")
    volume_ml: float | None = Field(None, description="Volume in mL")
    volume_l: float | None = Field(None, description="Volume in L")
    mw: float | None = Field(None, description="Molecular weight in g/mol")
    formula: str | None = Field(None, description="Molecular formula (e.g. 'NaCl', 'Ca(OH)2') — auto-calculates MW")


@mcp.tool()
async def calculate_molarity(params: CalcMolarityInput) -> str:
    """Calculate molarity/concentration. Provide 2-3 knowns (mass, moles, volume, MW, formula) → solves for unknowns. Handles M=n/V, n=m/MW, and all permutations."""
    try:
        result = calc_molarity(
            concentration=params.concentration_M,
            mass_grams=params.mass_grams,
            moles=params.moles,
            volume_ml=params.volume_ml,
            volume_l=params.volume_l,
            mw=params.mw,
            formula=params.formula,
        )
        lines = ["**Molarity Calculation**"]
        for k, v in result.items():
            label = k.replace("_", " ").replace("per", "/")
            lines.append(f"  {label}: {v}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


class CalcDilutionInput(BaseModel):
    """C₁V₁ = C₂V₂ solver. Provide 3 of 4 values → solves for the unknown."""
    c1: float | None = Field(None, description="Initial concentration")
    v1: float | None = Field(None, description="Initial volume")
    c2: float | None = Field(None, description="Final concentration")
    v2: float | None = Field(None, description="Final volume")
    c1_unit: str = Field("M", description="Unit for c1 (M, mM, μM, nM, mg/mL, %)")
    v1_unit: str = Field("mL", description="Unit for v1 (mL, L, μL)")
    c2_unit: str = Field("M", description="Unit for c2")
    v2_unit: str = Field("mL", description="Unit for v2")


@mcp.tool()
async def calculate_dilution(params: CalcDilutionInput) -> str:
    """C₁V₁ = C₂V₂ dilution calculator. Provide any 3 values → solves for the 4th. Reports solvent to add."""
    try:
        result = calc_dilution(
            c1=params.c1, v1=params.v1, c2=params.c2, v2=params.v2,
            c1_unit=params.c1_unit, v1_unit=params.v1_unit,
            c2_unit=params.c2_unit, v2_unit=params.v2_unit,
        )
        lines = ["**Dilution Calculation (C₁V₁ = C₂V₂)**"]
        for k, v in result.items():
            label = k.replace("_", " ")
            lines.append(f"  {label}: {v}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


class ReagentEntry(BaseModel):
    """A single reagent in a reaction mass calculation."""
    name: str = Field(description="Reagent name")
    mw: float | None = Field(None, description="Molecular weight (g/mol). Or provide formula.")
    formula: str | None = Field(None, description="Molecular formula — auto-calculates MW")
    equiv: float = Field(1.0, description="Equivalents (1.0 = stoichiometric)")
    mass_g: float | None = Field(None, description="Known mass in grams")
    moles: float | None = Field(None, description="Known moles")
    volume_ml: float | None = Field(None, description="Known volume in mL (requires density or molarity)")
    density: float | None = Field(None, description="Density in g/mL (for neat liquids)")
    molarity: float | None = Field(None, description="Solution concentration in M")
    limiting: bool = Field(False, description="Is this the limiting reagent?")


class CalcReactionMassInput(BaseModel):
    """Calculate masses/volumes for all reagents in a reaction from equivalents."""
    reagents: list[ReagentEntry] = Field(description="List of reagents with their properties")


@mcp.tool()
async def calculate_reaction_mass(params: CalcReactionMassInput) -> str:
    """Calculate mass, moles, and volume for each reagent in a reaction. Specify one limiting reagent with known mass/moles, set equivalents for others → get required amounts for all."""
    try:
        reagent_dicts = [r.model_dump(exclude_none=True) for r in params.reagents]
        result = calc_reaction_mass(reagent_dicts)
        lines = ["**Reaction Mass Calculation**"]
        if "limiting_reagent" in result:
            lines.append(f"Limiting reagent: {result['limiting_reagent']}")
            lines.append(f"Reference moles: {result.get('reference_moles', 'N/A')}")
        lines.append("")
        for entry in result.get("reagents", []):
            lines.append(f"**{entry.get('name', '?')}** ({entry.get('equiv', '?')} equiv)")
            for k in ["moles", "mass_g", "volume_ml"]:
                if k in entry:
                    lines.append(f"  {k.replace('_', ' ')}: {entry[k]:.4g}")
            lines.append("")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


class CalcYieldInput(BaseModel):
    """Calculate reaction yield."""
    actual_mass_g: float | None = Field(None, description="Actual product mass obtained (g)")
    actual_moles: float | None = Field(None, description="Actual product moles obtained")
    theoretical_mass_g: float | None = Field(None, description="Theoretical product mass (g)")
    theoretical_moles: float | None = Field(None, description="Theoretical product moles")
    product_mw: float | None = Field(None, description="Product molecular weight (g/mol)")
    product_formula: str | None = Field(None, description="Product formula — auto-calculates MW")
    limiting_mass_g: float | None = Field(None, description="Limiting reagent mass (g)")
    limiting_mw: float | None = Field(None, description="Limiting reagent MW (g/mol)")
    limiting_formula: str | None = Field(None, description="Limiting reagent formula")
    stoich_ratio: float = Field(1.0, description="Product:limiting reagent stoichiometric ratio")


@mcp.tool()
async def calculate_yield(params: CalcYieldInput) -> str:
    """Calculate percent yield. Provide actual product amount and either theoretical amount or limiting reagent info."""
    try:
        result = calc_yield(
            actual_mass_g=params.actual_mass_g,
            actual_moles=params.actual_moles,
            theoretical_mass_g=params.theoretical_mass_g,
            theoretical_moles=params.theoretical_moles,
            product_mw=params.product_mw,
            product_formula=params.product_formula,
            limiting_mass_g=params.limiting_mass_g,
            limiting_mw=params.limiting_mw,
            limiting_formula=params.limiting_formula,
            stoich_ratio=params.stoich_ratio,
        )
        lines = ["**Yield Calculation**"]
        for k, v in result.items():
            label = k.replace("_", " ").title()
            if isinstance(v, float):
                lines.append(f"  {label}: {v:.4g}")
            else:
                lines.append(f"  {label}: {v}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


class CalcConcentrationInput(BaseModel):
    """Convert between concentration units."""
    value: float = Field(description="Concentration value to convert")
    from_unit: str = Field(description="Source unit: M, mM, μM, nM, mg/mL, g/L, μg/mL, %w/v, %w/w, %v/v, ppm, ppb")
    to_unit: str = Field(description="Target unit (same options)")
    mw: float | None = Field(None, description="Molecular weight (required for molar ↔ mass conversions)")
    formula: str | None = Field(None, description="Molecular formula (auto-calculates MW)")
    density_solution: float = Field(1.0, description="Solution density in g/mL (for %w/w conversions)")


@mcp.tool()
async def calculate_concentration(params: CalcConcentrationInput) -> str:
    """Convert between concentration units (M, mM, μM, mg/mL, %w/v, ppm, ppb, etc.). MW required for molar↔mass conversions."""
    try:
        result = calc_concentration(
            value=params.value,
            from_unit=params.from_unit,
            to_unit=params.to_unit,
            mw=params.mw,
            formula=params.formula,
            density_solution=params.density_solution,
        )
        lines = [f"**Concentration Conversion**"]
        lines.append(f"  {params.value} {params.from_unit} = {result['result_value']:.6g} {params.to_unit}")
        if "intermediate_mg_per_mL" in result:
            lines.append(f"  (intermediate: {result['intermediate_mg_per_mL']:.4g} mg/mL)")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


# =============================================================================
# Bench Chemistry — Reference Lookups (7 tools, embedded knowledge)
# =============================================================================


class LookupReactionInput(BaseModel):
    """Search named organic chemistry reactions."""
    query: str = Field(description="Reaction name, category, substrate, or keyword (e.g. 'Suzuki', 'cross-coupling', 'aldehyde oxidation')")


@mcp.tool()
async def lookup_named_reaction(params: LookupReactionInput) -> str:
    """Search 86+ named organic reactions. Returns conditions, mechanism, substrate scope, limitations, and related reactions. Covers: cross-couplings, oxidations, reductions, rearrangements, cycloadditions, FGIs, radical, heterocycle synthesis."""
    results = lookup_reaction(params.query)
    if not results:
        return f"No named reactions found matching '{params.query}'. Try: Suzuki, Grignard, Wittig, Swern, Diels-Alder, aldol, etc."
    lines = []
    for r in results[:5]:
        lines.append(f"## {r['name']}")
        if r.get("aliases"):
            lines.append(f"**Also known as:** {', '.join(r['aliases'])}")
        lines.append(f"**Category:** {r['category']} — {r['type']}")
        lines.append(f"**Summary:** {r['summary']}")
        lines.append(f"**Conditions:** {r['conditions']}")
        lines.append(f"**Substrate scope:** {r['substrate']}")
        if r.get("limitations"):
            lines.append(f"**Limitations:** {r['limitations']}")
        if r.get("mechanism"):
            lines.append(f"**Mechanism:** {r['mechanism']}")
        if r.get("related"):
            lines.append(f"**Related:** {r['related']}")
        if r.get("refs"):
            lines.append(f"**Key reference:** {r['refs']}")
        lines.append("")
    if len(results) > 5:
        lines.append(f"*…and {len(results) - 5} more results*")
    return "\n".join(lines)


class LookupPGInput(BaseModel):
    """Search protecting groups."""
    query: str = Field(description="PG name or functional group (e.g. 'Boc', 'TBS', 'hydroxyl', 'amino')")
    functional_group: str | None = Field(None, description="Filter by: hydroxyl, amino, carbonyl, carboxyl")


@mcp.tool()
async def lookup_protecting_group(params: LookupPGInput) -> str:
    """Search 30 protecting groups across hydroxyl, amino, carbonyl, carboxyl. Returns protection/deprotection conditions and stability matrix (acid/base/nucleophile/oxidation/reduction/H₂-Pd). S=stable, M=moderate, L=labile."""
    results = lookup_pg(params.query, params.functional_group)
    if not results:
        return f"No protecting groups found for '{params.query}'. Try: Boc, Fmoc, TBS, PMB, MOM, Bn, etc."
    lines = []
    for pg in results[:10]:
        stab = pg.get("stability", {})
        stab_str = " | ".join(f"{k}:{v}" for k, v in stab.items())
        lines.append(f"## {pg['name']} ({pg.get('full_name', '')})")
        lines.append(f"**Protection:** {pg['protection']}")
        lines.append(f"**Deprotection:** {pg['deprotection']}")
        lines.append(f"**Stability:** {stab_str}")
        if pg.get("notes"):
            lines.append(f"**Notes:** {pg['notes']}")
        lines.append("")
    return "\n".join(lines)


class LookupWorkupInput(BaseModel):
    """Search workup procedures."""
    query: str = Field(description="Reaction type or workup keyword (e.g. 'LAH', 'Fieser', 'aqueous', 'extraction', 'Grignard')")


@mcp.tool()
async def lookup_workup_procedure(params: LookupWorkupInput) -> str:
    """Search bench workup procedures. Returns step-by-step protocols for: standard aqueous workup, acid-base extraction, Fieser (LAH/DIBAL), NH₄Cl quench, Rochelle's salt, and more."""
    results = lookup_workup(params.query)
    if not results:
        return f"No workup procedures found for '{params.query}'. Try: LAH, Fieser, extraction, aqueous, quench, Grignard."
    lines = []
    for w in results[:3]:
        lines.append(f"## {w['name']}")
        if w.get("use_when"):
            lines.append(f"**When to use:** {w['use_when']}")
        if w.get("steps"):
            lines.append("**Procedure:**")
            for i, step in enumerate(w["steps"], 1):
                lines.append(f"  {i}. {step}")
        if w.get("tips"):
            lines.append("**Tips:**")
            for tip in w["tips"]:
                lines.append(f"  • {tip}")
        lines.append("")
    return "\n".join(lines)


class LookupSolventInput(BaseModel):
    """Search solvent properties."""
    query: str = Field(description="Solvent name or property keyword (e.g. 'THF', 'DCM', 'polar aprotic', 'high boiling')")


@mcp.tool()
async def lookup_solvent_properties(params: LookupSolventInput) -> str:
    """Search 32 common lab solvents. Returns: boiling point, density, polarity index, dielectric constant, water miscibility, common uses, and safety notes."""
    results = lookup_solvent(params.query)
    if not results:
        return f"No solvents found for '{params.query}'. Try: THF, DCM, DMF, DMSO, EtOAc, hexanes, toluene, MeCN, etc."
    lines = []
    for s in results[:10]:
        lines.append(f"## {s['name']}")
        props = []
        if s.get("bp"):
            props.append(f"bp {s['bp']}°C")
        if s.get("density"):
            props.append(f"d={s['density']} g/mL")
        if s.get("polarity_index") is not None:
            props.append(f"polarity index={s['polarity_index']}")
        if s.get("dielectric") is not None:
            props.append(f"ε={s['dielectric']}")
        lines.append(f"**Properties:** {', '.join(props)}")
        if s.get("water_misc") is not None:
            lines.append(f"**Water miscible:** {'Yes' if s['water_misc'] else 'No'}")
        if s.get("uses"):
            lines.append(f"**Common uses:** {s['uses']}")
        if s.get("safety"):
            lines.append(f"**Safety:** {s['safety']}")
        lines.append("")
    return "\n".join(lines)


class LookupCoolingBathInput(BaseModel):
    """Find cooling bath recipe for a target temperature."""
    target_temp_c: float | None = Field(None, description="Target temperature in °C (e.g. -78, 0, -42)")


@mcp.tool()
async def lookup_cooling_bath(params: LookupCoolingBathInput) -> str:
    """Find cooling bath recipes. Provide target temp → get recipe. Covers -196°C (liq N₂) to +100°C. Common: -78°C dry ice/acetone, -42°C MeCN/dry ice, 0°C ice/water, -15°C ice/NaCl."""
    results = lookup_cooling_bath(params.target_temp_c)
    if not results:
        return "No cooling bath data available."
    lines = ["**Cooling Bath Recipes**"]
    for b in results:
        lines.append(f"  **{b['temp']}°C** — {b['recipe']}")
        if b.get("notes"):
            lines.append(f"    _{b['notes']}_")
    return "\n".join(lines)


class LookupTLCStainInput(BaseModel):
    """Search TLC stain recipes."""
    query: str = Field(description="Functional group or stain name (e.g. 'amine', 'aldehyde', 'KMnO4', 'CAM', 'ninhydrin', 'UV')")


@mcp.tool()
async def lookup_tlc_stain(params: LookupTLCStainInput) -> str:
    """Search TLC stain recipes by functional group or stain name. Returns: recipe, preparation, visualization procedure, and which functional groups each stain detects. Covers: KMnO₄, CAM, ninhydrin, anisaldehyde, vanillin, DNP, PMA, Dragendorff, I₂, UV."""
    results = lookup_tlc_stain(params.query)
    if not results:
        return f"No TLC stains found for '{params.query}'. Try: amine, alcohol, aldehyde, KMnO4, CAM, ninhydrin, UV, universal."
    lines = []
    for s in results[:8]:
        lines.append(f"## {s['name']}")
        if s.get("targets"):
            targets = s["targets"] if isinstance(s["targets"], str) else ", ".join(s["targets"])
            lines.append(f"**Detects:** {targets}")
        if s.get("recipe"):
            lines.append(f"**Recipe:** {s['recipe']}")
        if s.get("procedure"):
            lines.append(f"**Procedure:** {s['procedure']}")
        if s.get("color"):
            lines.append(f"**Color:** {s['color']}")
        lines.append("")
    return "\n".join(lines)


class LookupColumnInput(BaseModel):
    """Column chromatography guidance."""
    query: str = Field(description="Topic: 'solvent selection', 'loading', 'Rf', 'troubleshooting', 'streaking', 'gradient', 'dry loading', or specific problem")


@mcp.tool()
async def lookup_column_chromatography(params: LookupColumnInput) -> str:
    """Column chromatography guide. Search for: solvent selection, loading methods, Rf rules, troubleshooting (streaking, poor separation, compound stuck), gradient tips, fraction size, flow rate. Includes 8 common solvent systems."""
    results = lookup_column_guide(params.query)

    if isinstance(results, dict):
        lines = []
        for section, content in results.items():
            lines.append(f"## {section}")
            if isinstance(content, list):
                for item in content:
                    lines.append(f"  • {item}")
            elif isinstance(content, dict):
                for k, v in content.items():
                    lines.append(f"  **{k}:** {v}")
            else:
                lines.append(f"  {content}")
            lines.append("")
        # Append solvent systems
        lines.append("## Common Solvent Systems (Normal Phase)")
        for sys in CHROM_SOLVENT_SYSTEMS:
            lines.append(f"  • {sys.get('system', '')}: {sys.get('polarity', '')} — {sys.get('use', '')}")
        return "\n".join(lines)

    elif isinstance(results, list):
        return "\n".join(str(r) for r in results)
    else:
        return str(results)


# =============================================================================
# Peptide Chemistry — p2smi tools (6 tools)
# =============================================================================


class PeptideToSmilesInput(BaseModel):
    """Convert peptide sequence to SMILES."""
    sequence: str = Field(description="Amino acid sequence in 1-letter codes (supports 450+ AAs including NCAAs from SwissSidechain)")
    cyclization: str = Field("", description="Cyclization type: '' (linear), 'SS' (disulfide), 'HT' (head-to-tail), 'SCNT', 'SCCT', 'SCSC', or manual pattern like 'SSXXXCXXXCX'")


@mcp.tool()
async def peptide_to_smiles(params: PeptideToSmilesInput) -> str:
    """Convert peptide sequence to SMILES. Supports 450+ amino acids (canonical + noncanonical from SwissSidechain), D-stereochemistry, and 5 cyclization types: disulfide, head-to-tail, sidechain-to-sidechain, sidechain-to-N-term, sidechain-to-C-term."""
    try:
        result = sequence_to_smiles(params.sequence, params.cyclization)
        lines = [f"**Peptide → SMILES**"]
        lines.append(f"Sequence: {result['sequence']}")
        lines.append(f"Type: {result['type']}")
        if result.get('applied_constraint'):
            lines.append(f"Constraint: {result['applied_constraint']}")
        lines.append(f"SMILES: {result['smiles']}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error converting sequence: {e}"


class PeptideCyclizationInput(BaseModel):
    """Check cyclization options for a peptide."""
    sequence: str = Field(description="Amino acid sequence")


@mcp.tool()
async def peptide_cyclization_options(params: PeptideCyclizationInput) -> str:
    """Check which cyclization types a peptide sequence supports. Analyzes residues for disulfide-capable (Cys), nucleophilic sidechain (Lys, Ser, etc.), and electrophilic sidechain (Asp, Glu, etc.) positions."""
    try:
        result = get_cyclization_options(params.sequence)
        lines = [f"**Cyclization Options for {result['sequence']}** (length {result['length']})"]
        if result["num_options"] == 0:
            lines.append("No cyclization options available for this sequence.")
        else:
            lines.append(f"Found {result['num_options']} options:")
            for opt in result["cyclization_options"]:
                lines.append(f"  • **{opt['type']}** — pattern: {opt['constraint_pattern']}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


class GeneratePeptidesInput(BaseModel):
    """Generate random peptide sequences."""
    num_sequences: int = Field(10, description="Number of peptides to generate (max 100)")
    min_length: int = Field(8, description="Minimum sequence length")
    max_length: int = Field(20, description="Maximum sequence length")
    noncanonical_percent: float = Field(0.0, description="Fraction of noncanonical amino acids (0-1)")
    dextro_percent: float = Field(0.0, description="Fraction of D-amino acids (0-1)")
    cyclization: str = Field("none", description="Cyclization: 'none', 'all', or comma-separated: 'SS,HT,SCSC,SCNT,SCCT'")


@mcp.tool()
async def generate_peptide_library(params: GeneratePeptidesInput) -> str:
    """Generate random peptide sequences with controlled properties. Supports noncanonical amino acids (100+ from SwissSidechain), D-stereochemistry, and 5 cyclization types. Useful for library design and computational peptide chemistry."""
    try:
        result = generate_peptides(
            num_sequences=params.num_sequences,
            min_length=params.min_length,
            max_length=params.max_length,
            noncanonical_percent=params.noncanonical_percent,
            dextro_percent=params.dextro_percent,
            cyclization=params.cyclization,
        )
        lines = [f"**Generated {result['num_generated']} peptides**"]
        p = result["parameters"]
        lines.append(f"Length: {p['min_length']}–{p['max_length']}, NCAA: {p['noncanonical_percent']:.0%}, D-AA: {p['dextro_percent']:.0%}, cyclization: {p['cyclization']}")
        lines.append("")
        for seq in result["sequences"]:
            lines.append(f"  {seq['id']}: {seq['sequence']} (len={seq['length']})")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


class PeptidePropertiesInput(BaseModel):
    """Compute molecular properties for a peptide SMILES."""
    smiles: str = Field(description="Peptide SMILES string (use peptide_to_smiles to convert from sequence first)")


@mcp.tool()
async def peptide_properties(params: PeptidePropertiesInput) -> str:
    """Compute peptide molecular properties from SMILES: molecular weight, formula, logP, TPSA, H-bond donors/acceptors, rotatable bonds, ring count, fraction Csp3, heavy atoms, formal charge, and Lipinski Rule-of-5 evaluation."""
    try:
        result = get_peptide_properties(params.smiles)
        lines = ["**Peptide Molecular Properties**"]
        for k, v in result.items():
            if k == "SMILES":
                continue
            lines.append(f"  {k}: {v}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


class CheckSynthesisInput(BaseModel):
    """Check peptide synthesis feasibility."""
    sequence: str = Field(description="Amino acid sequence (natural amino acids, 1-letter codes)")


@mcp.tool()
async def check_peptide_synthesis(params: CheckSynthesisInput) -> str:
    """Evaluate solid-phase peptide synthesis (SPPS) feasibility. Checks: forbidden motifs (consecutive Pro, DG/DP, N/Q at N-term), cysteine overload, C-terminal Pro/Cys, glycine runs, length limits, hydrophobicity, and charge distribution."""
    try:
        result = check_synthesis_feasibility(params.sequence)
        lines = [f"**SPPS Synthesis Check: {result['sequence']}** (length {result['length']})"]
        lines.append(f"Verdict: **{result['verdict']}**")
        if result["issues"]:
            lines.append(f"Issues ({result['num_issues']}):")
            for issue in result["issues"]:
                lines.append(f"  ⚠ {issue}")
        else:
            lines.append("✓ No synthesis issues detected — sequence passes all checks.")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


class ModifyPeptideInput(BaseModel):
    """Apply chemical modifications to peptide SMILES."""
    smiles: str = Field(description="Peptide SMILES string")
    n_methylation: bool = Field(False, description="Apply random N-methylation to amide bonds")
    pegylation: bool = Field(False, description="Apply random PEGylation")
    methylation_fraction: float = Field(0.3, description="Fraction of amide sites to N-methylate (0-1)")


@mcp.tool()
async def modify_peptide(params: ModifyPeptideInput) -> str:
    """Apply N-methylation and/or PEGylation to a peptide SMILES. N-methylation improves metabolic stability and membrane permeability. PEGylation increases solubility and plasma half-life. Returns modified SMILES with validation."""
    try:
        result = modify_peptide_smiles(
            params.smiles,
            n_methylation=params.n_methylation,
            pegylation=params.pegylation,
            methylation_fraction=params.methylation_fraction,
        )
        lines = ["**Peptide Modification**"]
        lines.append(f"Modifications: {', '.join(result['modifications_applied']) or 'None'}")
        lines.append(f"Valid SMILES: {result['valid_smiles']}")
        if result["modified_smiles"]:
            lines.append(f"Modified SMILES: {result['modified_smiles']}")
        if result.get("error"):
            lines.append(f"Error: {result['error']}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}"


# =============================================================================
# Peptide Chemistry — pichemist pI calculation (1 tool)
# =============================================================================


class CalcPeptidePIInput(BaseModel):
    """Calculate peptide isoelectric point."""
    sequence: str | None = Field(None, description="Peptide sequence in 1-letter code (for canonical peptides)")
    smiles: str | None = Field(None, description="SMILES string (for modified/noncanonical peptides — pichemist cuts amide bonds and predicts pKas)")


@mcp.tool()
async def calculate_peptide_pi(params: CalcPeptidePIInput) -> str:
    """Calculate peptide isoelectric point (pI) using AstraZeneca's pichemist. Uses 8 different pKa reference sets (IPC2, ProMoST, Gauci, Grimsley, Thurlkill, Lehninger, Toseland) for consensus pI with error bars. Also returns charge at pH 7.4. Accepts FASTA sequence or SMILES (for modified peptides with noncanonical AAs)."""
    try:
        if params.sequence:
            result = calculate_pi_from_sequence(params.sequence)
        elif params.smiles:
            result = calculate_pi_from_smiles(params.smiles)
        else:
            return "Error: provide either sequence or smiles"

        if "error" in result:
            return f"Error: {result['error']}"

        lines = ["**Isoelectric Point (pI) Calculation** — pichemist (AstraZeneca)"]
        if result.get("sequence"):
            lines.append(f"Sequence: {result['sequence']}")
        if result.get("smiles"):
            lines.append(f"SMILES: {result['smiles'][:80]}...")

        lines.append(f"\n**Consensus pI: {result['pI_mean']:.2f} ± {result.get('pI_std', 0):.2f}**")
        lines.append(f"pI interval: {result.get('pI_interval', 'N/A')}")
        lines.append(f"Charge at pH 7.4: {result.get('charge_at_pH7_mean', 'N/A'):.2f}")

        pi_methods = result.get("pI_by_method", {})
        if pi_methods:
            lines.append("\npI by reference set:")
            for method, val in pi_methods.items():
                lines.append(f"  {method}: {val}")

        return "\n".join(lines)
    except Exception as e:
        return f"Error calculating pI: {e}"


# =============================================================================
# Peptide Chemistry — pep-calc.com API (1 tool, combines properties + MS)
# =============================================================================


class PepCalcInput(BaseModel):
    """Get peptide properties from pep-calc.com API."""
    sequence: str = Field(description="Peptide sequence (1-letter code, max 150 residues). Non-standard AAs in parentheses: (pS) for phosphoserine")
    n_term: str = Field("H", description="N-terminus: 'H' (free amine), 'Ac' (acetylated), etc.")
    c_term: str = Field("OH", description="C-terminus: 'OH' (free acid), 'NH2' (amidated), etc.")


@mcp.tool()
async def calculate_peptide_extinction(params: PepCalcInput) -> str:
    """Calculate peptide properties via pep-calc.com: MW, formula, pI, charge summary, and extinction coefficient at 280nm (oxidized/reduced). Supports 120+ amino acids including phospho-residues. Free API, no key required."""
    try:
        result = await pep_calc_properties(params.sequence, params.n_term, params.c_term)
        if result.get("basic_error"):
            return f"pep-calc.com API error: {result['basic_error']}. The API may be unreachable from this network. Use calculate_peptide_pi for pI and peptide_properties for MW/properties as alternatives."

        lines = ["**Peptide Properties** — pep-calc.com"]
        lines.append(f"Sequence: {params.sequence}")
        lines.append(f"Termini: N={params.n_term}, C={params.c_term}")
        for k, v in result.items():
            if v is not None and not k.endswith("_error"):
                label = k.replace("_", " ").title()
                lines.append(f"  {label}: {v}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}. The pep-calc.com API may be unreachable. Use peptide_properties and calculate_peptide_pi as local alternatives."


# =============================================================================
# Peptide MS — pep-calc.com ion series & peak assignment (2 tools)
# =============================================================================


class PepCalcIonsInput(BaseModel):
    """Get peptide fragment ion series for MS/MS."""
    sequence: str = Field(description="Peptide sequence (1-letter codes, max 150 residues)")
    n_term: str = Field("H", description="N-terminus: 'H' or 'Ac'")
    c_term: str = Field("OH", description="C-terminus: 'OH' or 'NH2'")


@mcp.tool()
async def get_peptide_ion_series(params: PepCalcIonsInput) -> str:
    """Get peptide MS/MS fragment ion series (b, y, a, c, z ions) for tandem mass spectrometry interpretation. Returns theoretical m/z values for each ion type at every cleavage position. Essential for de novo sequencing and spectral annotation."""
    try:
        result = await pep_calc_ion_series(params.sequence, params.n_term, params.c_term)
        if not result:
            return "pep-calc.com API returned no data. The API may be unreachable."
        lines = [f"**Peptide Ion Series** — {params.sequence}"]
        lines.append(f"Termini: N={params.n_term}, C={params.c_term}")
        for key, value in result.items():
            if isinstance(value, list) and value:
                lines.append(f"\n**{key} ions:**")
                for i, mz in enumerate(value, 1):
                    lines.append(f"  {key}{i}: {mz}")
            elif isinstance(value, (int, float)):
                lines.append(f"{key}: {value}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}. The pep-calc.com API may be unreachable."


class PepCalcMSAssignInput(BaseModel):
    """Assign MS peaks to peptide fragments."""
    sequence: str = Field(description="Peptide sequence (1-letter codes)")
    mz_values: list[float] = Field(description="Observed m/z values to assign")
    n_term: str = Field("H", description="N-terminus")
    c_term: str = Field("OH", description="C-terminus")


@mcp.tool()
async def assign_peptide_ms_peaks(params: PepCalcMSAssignInput) -> str:
    """Assign observed m/z peaks to peptide deletions, modifications, or fragment ions using pep-calc.com. Helps identify impurities, truncations, and post-translational modifications in MALDI/ESI-MS data."""
    try:
        result = await pep_calc_ms_assign(
            params.sequence, params.mz_values, params.n_term, params.c_term
        )
        if not result:
            return "pep-calc.com API returned no data."
        lines = [f"**MS Peak Assignment** — {params.sequence}"]
        lines.append(f"Query peaks: {params.mz_values}")
        if isinstance(result, dict):
            for key, value in result.items():
                if isinstance(value, list):
                    lines.append(f"\n**{key}:**")
                    for item in value[:20]:
                        lines.append(f"  {item}")
                else:
                    lines.append(f"{key}: {value}")
        elif isinstance(result, list):
            for item in result[:30]:
                lines.append(f"  {item}")
        return "\n".join(lines)
    except Exception as e:
        return f"Error: {e}. The pep-calc.com API may be unreachable."


# =============================================================================
# Bench Chemistry Expansion — buffers, amino acids, NMR solvents (3 tools)
# =============================================================================


class LookupBufferInput(BaseModel):
    """Search buffer recipes."""
    query: str = Field(description="Buffer name, pH value, category, or keyword (e.g. 'PBS', 'Tris', '7.4', 'electrophoresis', 'lysis', 'Western')")


@mcp.tool()
async def lookup_buffer_recipe(params: LookupBufferInput) -> str:
    """Search 20+ common laboratory buffer recipes. Returns: exact recipe with masses, preparation steps, pH range, pKa, and practical notes. Covers: PBS, TBS, Tris, HEPES, MES, MOPS, PIPES, citrate, acetate, carbonate, TAE, TBE, SDS-PAGE, transfer buffer, TE, RIPA lysis."""
    results = _lookup_buffer(params.query)
    if not results:
        return f"No buffer recipes found for '{params.query}'. Try: PBS, Tris, HEPES, TAE, TBE, SDS-PAGE, lysis, or a pH value."
    lines = []
    for buf in results[:5]:
        lines.append(f"## {buf['name']} (pH {buf['pH']})")
        lines.append(f"**Recipe:** {buf['recipe']}")
        lines.append(f"**Preparation:** {buf['preparation']}")
        lines.append(f"**pKa:** {buf.get('pKa', 'N/A')}, effective range: {buf.get('range', 'N/A')}")
        if buf.get("notes"):
            lines.append(f"**Notes:** {buf['notes']}")
        lines.append("")
    if len(results) > 5:
        lines.append(f"*…and {len(results) - 5} more results*")
    return "\n".join(lines)


class LookupAAInput(BaseModel):
    """Look up amino acid properties."""
    query: str = Field(description="1-letter code (e.g. 'C'), 3-letter code ('Cys'), full name ('Cysteine'), property class ('aromatic', 'polar', 'charged'), or 'all' for the full table")


@mcp.tool()
async def lookup_amino_acid_properties(params: LookupAAInput) -> str:
    """Look up amino acid properties: MW (residue), pKa values (α-carboxyl, α-amino, sidechain), pI, Kyte-Doolittle hydropathy index, structural class, and practical notes (UV absorption, phosphorylation, reactivity). Covers all 20 canonical amino acids."""
    results = _lookup_amino_acid(params.query)
    if not results:
        return f"No amino acids found for '{params.query}'. Use 1-letter code (C), 3-letter (Cys), name (Cysteine), class (aromatic), or 'all'."
    lines = []
    for aa in results[:20]:
        pka_side = f"{aa['pKa_side']}" if aa['pKa_side'] else "—"
        lines.append(f"**{aa['code1']} ({aa['code3']}) — {aa['name']}** [{aa['class']}]")
        lines.append(f"  MW(residue): {aa['mw']:.2f} Da | pI: {aa['pI']} | Hydropathy: {aa['hydropathy']}")
        lines.append(f"  pKₐ: α-COOH={aa['pKa_carboxyl']}, α-NH₂={aa['pKa_amino']}, sidechain={pka_side}")
        if aa.get("notes"):
            lines.append(f"  _{aa['notes']}_")
        lines.append("")
    return "\n".join(lines)


class LookupNMRSolventInput(BaseModel):
    """Search NMR solvent reference."""
    query: str = Field(description="Solvent name or abbreviation (e.g. 'CDCl3', 'DMSO', 'D2O', 'methanol') or 'all' for complete table")


@mcp.tool()
async def lookup_nmr_solvent(params: LookupNMRSolventInput) -> str:
    """NMR solvent reference: residual ¹H and ¹³C chemical shifts, water peak position, multiplicity, boiling point, density. Covers: CDCl₃, DMSO-d₆, MeOD, D₂O, acetone-d₆, MeCN-d₃, C₆D₆, THF-d₈, CD₂Cl₂, DMF-d₇, pyridine-d₅, TFA-d. Essential for spectral interpretation and thesis writing."""
    results = _lookup_nmr_solvent(params.query)
    if not results:
        return f"No NMR solvents found for '{params.query}'. Try: CDCl3, DMSO, D2O, MeOD, acetone, benzene, THF, or 'all'."
    lines = []
    for sol in results[:12]:
        lines.append(f"## {sol['name']} ({sol['abbrev']})")
        lines.append(f"  ¹H residual: {sol['h_residual']} ppm")
        if sol.get("c_residual"):
            lines.append(f"  ¹³C residual: {sol['c_residual']} ppm ({sol.get('c_multiplicity', '')})")
        if sol.get("water_peak"):
            lines.append(f"  Water peak: {sol['water_peak']} ppm")
        lines.append(f"  bp: {sol['bp']}°C, d={sol['density']} g/mL")
        if sol.get("notes"):
            lines.append(f"  _{sol['notes']}_")
        lines.append("")
    return "\n".join(lines)


# =============================================================================
# Chemistry Utilities — isotope pattern, CAS, units, periodic table, pH (6 tools)
# =============================================================================


class IsotopePatternInput(BaseModel):
    """Calculate isotope distribution."""
    formula: str = Field("", description="Molecular formula (e.g. 'C9H8O4', 'C14H19N3O4S')")
    smiles: str = Field("", description="SMILES string (alternative — formula auto-extracted)")
    charge: int = Field(1, description="Charge state for m/z calculation (default 1 for [M+H]⁺)")
    top_n: int = Field(8, description="Number of peaks to return")


@mcp.tool()
async def calculate_isotope_pattern(params: IsotopePatternInput) -> str:
    """Calculate theoretical isotope distribution from molecular formula or SMILES. Returns monoisotopic mass, average mass, and isotope envelope (m/z + relative abundance). Essential for mass spec data interpretation, especially for compounds containing Cl, Br, S."""
    result = _calc_isotope_pattern(
        formula=params.formula, smiles=params.smiles,
        charge=params.charge, top_n=params.top_n,
    )
    if result.get("error"):
        return f"Error: {result['error']}"
    lines = [f"**Isotope Pattern — {result['formula']}** (z={result['charge']})"]
    lines.append(f"Monoisotopic mass: {result['monoisotopic_mass']} Da")
    lines.append(f"Monoisotopic m/z: {result['monoisotopic_mz']}")
    lines.append(f"Average mass: {result['average_mass']} Da")
    lines.append("\nIsotope envelope:")
    for p in result["pattern"]:
        bar = "█" * int(p["relative"] / 2)
        lines.append(f"  m/z {p['mz']:.4f}: {p['relative']:6.2f}% {bar}")
    return "\n".join(lines)


class ValidateCASInput(BaseModel):
    """Validate CAS registry number."""
    cas_number: str = Field(description="CAS number (e.g. '50-78-2' for aspirin, or '50782')")


@mcp.tool()
async def validate_cas_number(params: ValidateCASInput) -> str:
    """Validate a CAS registry number by checking the check digit. Returns whether the number is valid and the properly formatted CAS string. Useful for verifying catalog numbers and safety data lookup."""
    result = _validate_cas(params.cas_number)
    if result["valid"]:
        return f"✓ **{result['cas']}** is a valid CAS registry number."
    else:
        return f"✗ **{result['cas']}** is NOT valid: {result['error']}"


class ConvertUnitsInput(BaseModel):
    """Convert between scientific units."""
    value: float = Field(description="Numeric value to convert")
    from_unit: str = Field(description="Source unit (e.g. 'kcal/mol', 'Å', '°C', 'mg', 'mL', 'atm', 'mmol')")
    to_unit: str = Field(description="Target unit (must be same category as source)")


@mcp.tool()
async def convert_units(params: ConvertUnitsInput) -> str:
    """Convert between scientific units. Supports: mass (g↔mg↔μg↔Da↔kDa), volume (L↔mL↔μL), length (m↔nm↔Å↔pm), energy (J↔kJ↔kcal↔kcal/mol↔kJ/mol↔eV↔hartree↔cm⁻¹), pressure (Pa↔atm↔bar↔torr↔mmHg), time (s↔ms↔μs↔ns↔min↔h), amount (mol↔mmol↔μmol), temperature (°C↔°F↔K)."""
    result = _convert_units(params.value, params.from_unit, params.to_unit)
    if result.get("error"):
        return f"Error: {result['error']}"
    r = result["result"]
    # Smart formatting
    if abs(r) >= 1000 or (abs(r) < 0.001 and r != 0):
        r_str = f"{r:.6e}"
    else:
        r_str = f"{r:.6g}"
    return f"**{params.value} {params.from_unit}** = **{r_str} {params.to_unit}** ({result['category']})"


class LookupElementInput(BaseModel):
    """Look up periodic table element."""
    query: str = Field(description="Element symbol (Fe), name (Iron), or atomic number (26)")


@mcp.tool()
async def lookup_periodic_table(params: LookupElementInput) -> str:
    """Look up element properties from the periodic table: atomic number, atomic mass, electron configuration, electronegativity (Pauling), group, period, block, and category. Covers all commonly encountered elements in organic, inorganic, and materials chemistry."""
    result = _lookup_element(params.query)
    if result is None:
        return f"Element '{params.query}' not found. Try symbol (Fe), name (Iron), or atomic number (26)."
    if isinstance(result, list):
        lines = [f"Multiple matches for '{params.query}':"]
        for el in result[:10]:
            lines.append(f"  **{el['symbol']}** (Z={el['Z']}) — {el['name']}: {el['mass']} g/mol")
        return "\n".join(lines)
    el = result
    lines = [f"## {el['symbol']} — {el['name']}"]
    lines.append(f"  Atomic number: {el['Z']}")
    lines.append(f"  Atomic mass: {el['mass']} g/mol")
    lines.append(f"  Electron config: {el['config']}")
    en = f"{el['en']}" if el['en'] else "N/A"
    lines.append(f"  Electronegativity (Pauling): {en}")
    lines.append(f"  Group: {el.get('group', 'N/A')}, Period: {el['period']}, Block: {el['block']}")
    lines.append(f"  Category: {el['category']}")
    return "\n".join(lines)


class CalcBufferPHInput(BaseModel):
    """Henderson-Hasselbalch buffer pH calculator."""
    pKa: float | None = Field(None, description="pKa of the buffer acid/base pair")
    buffer_name: str | None = Field(None, description="Buffer name to look up pKa (e.g. 'Tris', 'HEPES', 'phosphate_2', 'acetate')")
    acid_conc: float | None = Field(None, description="Concentration of acid form [HA]")
    base_conc: float | None = Field(None, description="Concentration of conjugate base form [A⁻]")
    ratio_base_acid: float | None = Field(None, description="[A⁻]/[HA] ratio (alternative to concentrations)")
    target_ph: float | None = Field(None, description="Target pH → calculates required ratio")


@mcp.tool()
async def calculate_buffer_ph(params: CalcBufferPHInput) -> str:
    """Henderson-Hasselbalch calculator. Two modes: (1) Given pKa + concentrations → calculate pH. (2) Given pKa + target pH → calculate required acid:base ratio. Knows pKa for 20+ common buffers: Tris (8.07), HEPES (7.55), phosphate (7.20), MOPS (7.20), MES (6.15), acetate (4.76), etc."""
    result = _calc_buffer_ph(
        pKa=params.pKa,
        buffer_name=params.buffer_name,
        acid_conc=params.acid_conc,
        base_conc=params.base_conc,
        ratio_base_acid=params.ratio_base_acid,
        target_ph=params.target_ph,
    )
    if result.get("error"):
        return f"Error: {result['error']}"
    lines = [f"**Buffer pH Calculation** (Henderson-Hasselbalch)"]
    lines.append(f"pKa: {result['pKa']}")
    if result.get("buffer"):
        lines.append(f"Buffer: {result['buffer']}")
    if result.get("species"):
        lines.append(f"Species: {result['species']}")
    lines.append(f"Effective range: {result['buffer_range']}")
    if "calculated_pH" in result:
        lines.append(f"\n**Calculated pH: {result['calculated_pH']}**")
        lines.append(f"[A⁻]/[HA] ratio: {result['ratio_base_acid']}")
    if "target_pH" in result:
        lines.append(f"\nTarget pH: {result['target_pH']}")
        lines.append(f"**Required [A⁻]/[HA] ratio: {result['required_ratio_base_acid']}**")
        if result.get("note"):
            lines.append(result["note"])
    if result.get("warning"):
        lines.append(f"\n⚠️ {result['warning']}")
    return "\n".join(lines)


# =============================================================================
# Writing & Publication Tools (10 tools)
# =============================================================================

# Tool 69: format_citation — DOI → formatted reference
# =============================================================================


class FormatCitationInput(BaseModel):
    """Format a citation from DOI."""
    doi: str = Field(description="DOI to format (e.g. '10.1038/s41586-020-2649-2')")
    style: str = Field(
        default="acs",
        description="Citation style: 'acs' (ACS), 'apa' (APA 7th), 'nature', 'vancouver', 'chicago', 'mla', 'bibtex'"
    )


@mcp.tool()
async def format_citation(params: FormatCitationInput) -> str:
    """Format a DOI into a properly styled citation. Styles: ACS (default), APA 7th, Nature, Vancouver, Chicago, MLA, BibTeX. Uses Crossref content negotiation for authoritative metadata. Essential for manuscripts, theses, and reference lists."""
    result = await _format_citation(params.doi, style=params.style)
    if result.get("error"):
        return f"Error: {result['error']}"
    lines = [f"**{params.style.upper()} format:**\n"]
    lines.append(result.get("formatted", ""))
    if result.get("doi"):
        lines.append(f"\nDOI: https://doi.org/{result['doi']}")
    return "\n".join(lines)


# =============================================================================
# Tool 70: build_bibliography — batch DOIs → formatted reference list
# =============================================================================


class BuildBibliographyInput(BaseModel):
    """Build a formatted bibliography."""
    dois: list[str] = Field(description="List of DOIs to format")
    style: str = Field(
        default="acs",
        description="Citation style: 'acs', 'apa', 'nature', 'vancouver', 'chicago', 'mla', 'bibtex'"
    )
    numbered: bool = Field(default=True, description="Number the references")


@mcp.tool()
async def build_bibliography(params: BuildBibliographyInput) -> str:
    """Build a formatted bibliography from a list of DOIs. Outputs a numbered reference list in the specified style (ACS, APA, Nature, etc.). Ideal for thesis chapters, manuscript reference sections, and literature reviews."""
    result = await _build_bibliography(params.dois, style=params.style)
    if result.get("error"):
        return f"Error: {result['error']}"
    entries = result.get("entries", [])
    if not entries:
        return "No citations could be formatted."
    lines = [f"**Bibliography ({params.style.upper()}, {len(entries)} references)**\n"]
    for i, entry in enumerate(entries, 1):
        prefix = f"{i}. " if params.numbered else "• "
        if entry.get("error"):
            lines.append(f"{prefix}[Error for {entry.get('doi', '?')}]: {entry['error']}")
        else:
            lines.append(f"{prefix}{entry.get('formatted', '')}")
    if result.get("errors"):
        lines.append(f"\n⚠️ {result['errors']} DOI(s) could not be resolved.")
    return "\n".join(lines)


# =============================================================================
# Tool 71: lookup_iupac_name — SMILES → IUPAC name
# =============================================================================


class IUPACNameInput(BaseModel):
    """Look up IUPAC name from SMILES."""
    smiles: str = Field(description="SMILES string (e.g. 'CC(=O)Oc1ccccc1C(=O)O' for aspirin)")


@mcp.tool()
async def lookup_iupac_name(params: IUPACNameInput) -> str:
    """Convert a SMILES string to its IUPAC systematic name using PubChem. Returns the preferred IUPAC name along with any common synonyms. Useful for writing experimental sections and compound characterization."""
    result = await _iupac_from_smiles(params.smiles)
    if result.get("error"):
        return f"Error: {result['error']}"
    lines = []
    if result.get("iupac_name"):
        lines.append(f"**IUPAC Name:** {result['iupac_name']}")
    if result.get("cid"):
        lines.append(f"PubChem CID: {result['cid']}")
    if result.get("molecular_formula"):
        lines.append(f"Formula: {result['molecular_formula']}")
    if result.get("synonyms"):
        syns = result["synonyms"][:5]
        lines.append(f"Synonyms: {', '.join(syns)}")
    return "\n".join(lines) if lines else "No IUPAC name found for this SMILES."


# =============================================================================
# Tool 72: name_to_smiles — compound name → SMILES
# =============================================================================


class NameToSMILESInput(BaseModel):
    """Convert compound name to SMILES."""
    name: str = Field(description="Compound name (e.g. 'caffeine', 'ibuprofen', 'Fmoc-Gly-OH')")


@mcp.tool()
async def name_to_smiles(params: NameToSMILESInput) -> str:
    """Convert a compound name to its SMILES string using PubChem. Returns canonical and isomeric SMILES, InChI, InChIKey, molecular formula, and MW. Useful for quickly getting machine-readable structures from common names, trade names, or reagent labels."""
    result = await _smiles_from_name(params.name)
    if result.get("error"):
        return f"Error: {result['error']}"
    lines = [f"**{params.name}**\n"]
    if result.get("canonical_smiles"):
        lines.append(f"Canonical SMILES: `{result['canonical_smiles']}`")
    if result.get("isomeric_smiles"):
        lines.append(f"Isomeric SMILES: `{result['isomeric_smiles']}`")
    if result.get("inchi"):
        lines.append(f"InChI: `{result['inchi']}`")
    if result.get("inchikey"):
        lines.append(f"InChIKey: `{result['inchikey']}`")
    if result.get("molecular_formula"):
        lines.append(f"Formula: {result['molecular_formula']}")
    if result.get("molecular_weight"):
        lines.append(f"MW: {result['molecular_weight']} g/mol")
    if result.get("cid"):
        lines.append(f"PubChem CID: {result['cid']}")
    return "\n".join(lines)


# =============================================================================
# Tool 73: format_molecular_formula — formula → LaTeX/HTML/Unicode
# =============================================================================


class FormatFormulaInput(BaseModel):
    """Format molecular formula with proper subscripts."""
    formula: str = Field(description="Molecular formula (e.g. 'C6H12O6', 'CH3CO2H', 'Fe2O3')")
    output_format: str = Field(
        default="all",
        description="Output format: 'unicode' (C₆H₁₂O₆), 'latex' (\\ce{...}), 'html' (<sub>), or 'all'"
    )


@mcp.tool()
async def format_molecular_formula(params: FormatFormulaInput) -> str:
    """Format a molecular formula with proper subscript notation. Outputs in Unicode (C₆H₁₂O₆), LaTeX (\\ce{C6H12O6}), HTML (<sub>6</sub>), or all formats. Handles complex formulae including parentheses and charges. Essential for manuscripts, theses, and presentations."""
    if params.output_format == "all":
        results = {}
        for fmt in ("unicode", "latex", "html"):
            r = _format_molecular_formula(params.formula, output_format=fmt)
            results[fmt] = r
        lines = [f"**{params.formula}** — formatted:\n"]
        u = results.get("unicode", {})
        if u.get("unicode"):
            lines.append(f"Unicode:  {u['unicode']}")
        la = results.get("latex", {})
        if la.get("latex_inline"):
            lines.append(f"LaTeX:    `{la['latex_inline']}`")
        if la.get("latex_ce"):
            lines.append(f"mhchem:   `{la['latex_ce']}`")
        h = results.get("html", {})
        if h.get("html"):
            lines.append(f"HTML:     `{h['html']}`")
        return "\n".join(lines)
    else:
        result = _format_molecular_formula(params.formula, output_format=params.output_format)
        if result.get("error"):
            return f"Error: {result['error']}"
        lines = [f"**{params.formula}** ({params.output_format}):\n"]
        for key, val in result.items():
            if key != "formula" and key != "format":
                lines.append(f"{key}: `{val}`" if "{" in str(val) or "<" in str(val) else f"{key}: {val}")
        return "\n".join(lines)


# =============================================================================
# Tool 74: lookup_experimental_template — reaction type → experimental section
# =============================================================================


class ExperimentalTemplateInput(BaseModel):
    """Look up experimental section template."""
    query: str = Field(
        description="Reaction type (e.g. 'suzuki', 'amide coupling', 'grignard', 'hydrogenation', "
        "'reductive amination', 'esterification', 'wittig', 'peptide coupling') or 'all' for the full list"
    )


@mcp.tool()
async def lookup_experimental_template(params: ExperimentalTemplateInput) -> str:
    """Get a template for writing the Experimental section of a paper or thesis. Covers 18 common reaction types (cross-coupling, oxidation, reduction, heterocycle synthesis, etc.) with fill-in-the-blank templates that follow journal conventions. Each template includes: reaction setup, workup, purification, and characterization data format."""
    result = _lookup_experimental_template(params.query)
    if result is None:
        return f"No template found for '{params.query}'. Try 'all' to see available templates."
    if isinstance(result, list):
        # Multiple matches or listing
        lines = [f"**Experimental templates matching '{params.query}'** ({len(result)} results):\n"]
        for tmpl in result:
            lines.append(f"- **{tmpl['name']}** ({tmpl.get('category', 'General')})")
        lines.append("\nRequest a specific template by name for the full text.")
        return "\n".join(lines)
    # Single result
    t = result
    lines = [f"## {t['name']}"]
    if t.get("category"):
        lines.append(f"Category: {t['category']}\n")
    lines.append("### Template\n")
    lines.append(t.get("template", ""))
    if t.get("notes"):
        lines.append(f"\n### Notes\n{t['notes']}")
    if t.get("safety"):
        lines.append(f"\n⚠️ **Safety:** {t['safety']}")
    return "\n".join(lines)


# =============================================================================
# Tool 75: lookup_journal_guide — journal → submission formatting requirements
# =============================================================================


class JournalGuideInput(BaseModel):
    """Look up journal submission guide."""
    query: str = Field(
        description="Journal name or abbreviation (e.g. 'JACS', 'Angew. Chem.', 'Nature Chemistry', "
        "'JOC', 'Org. Lett.') or 'all' for the full list"
    )


@mcp.tool()
async def lookup_journal_guide(params: JournalGuideInput) -> str:
    """Look up submission formatting guidelines for major chemistry journals. Covers 12 top journals (JACS, Angew. Chem., Nature Chemistry, JOC, Org. Lett., Chemical Science, etc.). Returns: citation style, word/page limits, figure requirements, SI guidelines, common mistakes, and submission tips."""
    result = _lookup_journal_guide(params.query)
    if result is None:
        return f"No guide found for '{params.query}'. Try 'all' to see available journals."
    if isinstance(result, list):
        lines = [f"**Journal guides matching '{params.query}'** ({len(result)} results):\n"]
        for g in result:
            name = g.get("name", "?")
            pub = g.get("publisher", "")
            lines.append(f"- **{name}** ({pub})")
        lines.append("\nRequest a specific journal for the full guide.")
        return "\n".join(lines)
    g = result
    lines = [f"## {g['name']}"]
    if g.get("publisher"):
        lines.append(f"Publisher: {g['publisher']}")
    if g.get("issn"):
        lines.append(f"ISSN: {g['issn']}")
    lines.append("")
    if g.get("citation_style"):
        lines.append(f"**Citation style:** {g['citation_style']}")
    if g.get("word_limit"):
        lines.append(f"**Word limit:** {g['word_limit']}")
    if g.get("page_limit"):
        lines.append(f"**Page limit:** {g['page_limit']}")
    if g.get("abstract_limit"):
        lines.append(f"**Abstract limit:** {g['abstract_limit']}")
    if g.get("figure_guidelines"):
        lines.append(f"\n**Figures:** {g['figure_guidelines']}")
    if g.get("si_guidelines"):
        lines.append(f"**SI:** {g['si_guidelines']}")
    if g.get("submission_format"):
        lines.append(f"**Format:** {g['submission_format']}")
    if g.get("common_mistakes"):
        lines.append("\n**Common mistakes:**")
        for m in g["common_mistakes"]:
            lines.append(f"  • {m}")
    if g.get("tips"):
        lines.append("\n**Tips:**")
        for tip in g["tips"]:
            lines.append(f"  • {tip}")
    if g.get("url"):
        lines.append(f"\nAuthor guidelines: {g['url']}")
    return "\n".join(lines)


# =============================================================================
# Tool 76: generate_si_checklist — Supporting Information checklist
# =============================================================================


class SIChecklistInput(BaseModel):
    """Generate a Supporting Information checklist."""
    compound_type: str = Field(
        default="small molecule",
        description="Type: 'small molecule', 'peptide', 'polymer', 'material', 'natural product'"
    )
    content_types: list[str] | None = Field(
        default=None,
        description="Analytical data included (e.g. ['1h_nmr', '13c_nmr', 'hrms', 'hplc', 'ir', 'mp', 'optical_rotation', 'x-ray']). None → standard minimum."
    )


@mcp.tool()
async def generate_si_checklist(params: SIChecklistInput) -> str:
    """Generate a Supporting Information (SI) checklist for a chemistry publication. Tailored to compound type (small molecule, peptide, polymer, material, natural product). Lists required characterization data, formatting guidelines, and common reviewer complaints. Saves time during submission preparation."""
    result = _get_si_checklist(
        content_types=params.content_types,
        compound_type=params.compound_type,
    )
    if result.get("error"):
        return f"Error: {result['error']}"
    lines = [f"## SI Checklist — {result.get('compound_type', params.compound_type).title()}\n"]
    checklist = result.get("checklist", [])
    for item in checklist:
        if isinstance(item, dict):
            check = "☐"
            lines.append(f"{check} **{item.get('item', '')}**: {item.get('details', '')}")
        else:
            lines.append(f"☐ {item}")
    if result.get("general_tips"):
        lines.append("\n**General tips:**")
        for tip in result["general_tips"]:
            lines.append(f"  • {tip}")
    return "\n".join(lines)


# =============================================================================
# Tool 77: lookup_abbreviation — chemistry abbreviation reference
# =============================================================================


class LookupAbbreviationInput(BaseModel):
    """Look up chemistry abbreviations."""
    query: str = Field(
        description="Abbreviation (e.g. 'THF', 'DCM', 'TEA', 'DIPEA', 'HATU', 'EDC') "
        "or category ('solvents', 'reagents', 'protecting_groups', 'techniques', 'units', 'spectroscopy', 'general')"
    )


@mcp.tool()
async def lookup_abbreviation(params: LookupAbbreviationInput) -> str:
    """Look up chemistry abbreviations and acronyms. 193 entries across 7 categories: solvents, reagents, protecting groups, techniques, units, spectroscopy, and general chemistry. Returns full expansion. Also accepts category names to list all abbreviations in that category. Essential for decoding literature and ensuring correct abbreviation usage in manuscripts."""
    # First try as a specific abbreviation
    result = _lookup_abbreviation(params.query)
    if result.get("num_results", 0) > 0:
        lines = []
        for abbr, data in result.get("results", {}).items():
            lines.append(f"**{abbr}** = {data['meaning']}  ({data.get('category', '')})")
        return "\n".join(lines)
    # Try as category
    result = _get_abbreviations(params.query)
    if isinstance(result, dict) and not result.get("error"):
        lines = [f"**{params.query.title()} abbreviations:**\n"]
        for cat, abbrevs in result.items():
            if isinstance(abbrevs, dict):
                for abbr, meaning in sorted(abbrevs.items()):
                    lines.append(f"  **{abbr}** = {meaning}")
        return "\n".join(lines) if len(lines) > 1 else f"No abbreviations found for '{params.query}'."
    return f"Abbreviation '{params.query}' not found. Try a specific abbreviation or category: solvents, reagents, protecting_groups, techniques, units, spectroscopy, general."


# =============================================================================
# Tool 78: get_thesis_guide — thesis section writing guidance
# =============================================================================


class ThesisGuideInput(BaseModel):
    """Get thesis writing guidance."""
    section: str = Field(
        description="Thesis section: 'introduction', 'literature_review', 'methods', 'results', "
        "'discussion', 'conclusion', 'abstract', 'acknowledgements', or 'all' for overview"
    )


@mcp.tool()
async def get_thesis_guide(params: ThesisGuideInput) -> str:
    """Get writing guidance for each section of a chemistry thesis or dissertation. Covers: introduction, literature review, methods/experimental, results, discussion, conclusion, abstract, and acknowledgements. For each section: purpose, structure, writing tips, and common mistakes to avoid."""
    if params.section.lower() == "all":
        sections = ["abstract", "introduction", "literature_review", "methods",
                     "results", "discussion", "conclusion", "acknowledgements"]
        lines = ["## Thesis Writing Guide — Section Overview\n"]
        for s in sections:
            guide = _get_thesis_guide(s)
            if guide:
                lines.append(f"**{guide.get('name', s.replace('_', ' ').title())}** — {guide.get('purpose', '')}")
        lines.append("\nRequest a specific section for detailed guidance.")
        return "\n".join(lines)
    result = _get_thesis_guide(params.section)
    if result is None:
        return f"No guide found for '{params.section}'. Available: introduction, literature_review, methods, results, discussion, conclusion, abstract, acknowledgements."
    lines = [f"## {result.get('name', params.section.replace('_', ' ').title())}"]
    if result.get("purpose"):
        lines.append(f"\n**Purpose:** {result['purpose']}")
    if result.get("structure"):
        lines.append("\n**Structure:**")
        for s in result["structure"]:
            if isinstance(s, dict):
                lines.append(f"  {s.get('order', '')}. **{s.get('name', '')}** — {s.get('content', '')}")
            else:
                lines.append(f"  • {s}")
    if result.get("tips"):
        lines.append("\n**Tips:**")
        for tip in result["tips"]:
            lines.append(f"  • {tip}")
    if result.get("common_mistakes"):
        lines.append("\n**Common mistakes:**")
        for m in result["common_mistakes"]:
            lines.append(f"  ✗ {m}")
    return "\n".join(lines)


# =============================================================================
# Entry point
# =============================================================================


def main():
    """Run the MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
