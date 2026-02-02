"""
scholarly-mcp: Multi-source academic research MCP server.

Provides Claude with access to:
  - Literature search & paper metadata (Crossref, OpenAlex, Semantic Scholar)
  - Citation graphs & influence analysis (Semantic Scholar)
  - Paper recommendations (Semantic Scholar)
  - Author profiles & bibliometrics (OpenAlex, Semantic Scholar)
  - Research topic trend analysis (OpenAlex)
  - Open access PDF discovery (Unpaywall)
  - ChemRxiv preprint search (Crossref + OpenAlex)
  - Chemical compound lookup (PubChem, Common Chemistry/CAS)
  - Web of Science search (optional, credential-gated)
"""

from __future__ import annotations

import asyncio
import logging
from typing import Optional

from pydantic import BaseModel, Field
from mcp.server.fastmcp import FastMCP

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
# Entry point
# =============================================================================


def main():
    """Run the MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
