"""
API clients for scholarly-mcp.

All functions are async, return parsed dicts/lists or None on failure.
Each API module is self-contained with its own headers and error handling.
"""

from __future__ import annotations

import asyncio
import contextlib
import logging
import os
import re
from typing import Any

import httpx

logger = logging.getLogger("scholarly-mcp")

# =============================================================================
# Configuration (from environment variables)
# =============================================================================

VERSION = "3.0.0"
USER_AGENT = (
    f"scholarly-mcp/{VERSION} "
    "(https://github.com/JonasRackl/chemrxiv-mcp; "
    "mailto:scholarly-mcp@users.noreply.github.com)"
)
TIMEOUT = 30

# Optional credentials — features activate when set
S2_API_KEY: str | None = os.environ.get("S2_API_KEY")
OPENALEX_EMAIL: str | None = os.environ.get("OPENALEX_EMAIL")
UNPAYWALL_EMAIL: str | None = os.environ.get(
    "UNPAYWALL_EMAIL",
    os.environ.get("OPENALEX_EMAIL", "scholarly-mcp@users.noreply.github.com"),
)
WOS_API_KEY: str | None = os.environ.get("WOS_API_KEY")
MP_API_KEY: str | None = os.environ.get("MP_API_KEY")
RXN_API_KEY: str | None = os.environ.get("RXN_API_KEY")
COMPTOX_API_KEY: str | None = os.environ.get("COMPTOX_API_KEY")

# Base URLs
CROSSREF_BASE = "https://api.crossref.org"
OPENALEX_BASE = "https://api.openalex.org"
S2_BASE = "https://api.semanticscholar.org"
UNPAYWALL_BASE = "https://api.unpaywall.org/v2"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CAS_BASE = "https://commonchemistry.cas.org/api"
WOS_BASE = "https://wos-api.clarivate.com/api/wos"
NIST_BASE = "https://webbook.nist.gov/cgi/cbook.cgi"
MP_BASE = "https://api.materialsproject.org"
RXN_BASE = "https://rxn.res.ibm.com/rxn/api/api/v1"
UNICHEM_BASE = "https://www.ebi.ac.uk/unichem/api/v1"
COD_BASE = "https://www.crystallography.net/cod"
COMPTOX_BASE = "https://api-ccte.epa.gov"
MASSBANK_BASE = "https://massbank.eu/MassBank/api"
BINDINGDB_BASE = "https://bindingdb.org/axis2/services/BDBService"
PDB_DATA_BASE = "https://data.rcsb.org/rest/v1/core"
PDB_SEARCH_BASE = "https://search.rcsb.org/rcsbsearch/v2/query"
NPCLASSIFIER_BASE = "https://npclassifier.gnps2.org"

# Semantic Scholar field sets
S2_SEARCH_FIELDS = (
    "paperId,externalIds,title,abstract,year,venue,citationCount,"
    "influentialCitationCount,isOpenAccess,openAccessPdf,tldr,authors"
)
S2_DETAIL_FIELDS = (
    "paperId,externalIds,url,title,abstract,venue,publicationVenue,year,"
    "referenceCount,citationCount,influentialCitationCount,isOpenAccess,"
    "openAccessPdf,fieldsOfStudy,s2FieldsOfStudy,publicationTypes,"
    "publicationDate,journal,authors,tldr"
)
S2_CITATION_FIELDS = (
    "paperId,title,year,citationCount,authors,intents,isInfluential,contexts"
)
S2_AUTHOR_FIELDS = (
    "authorId,name,affiliations,paperCount,citationCount,hIndex"
)
S2_AUTHOR_DETAIL_FIELDS = (
    "authorId,name,affiliations,paperCount,citationCount,hIndex,"
    "papers,papers.paperId,papers.title,papers.year,papers.citationCount,"
    "papers.venue,papers.externalIds"
)

# ChemRxiv DOI prefix (for ChemRxiv-specific searches via Crossref)
CHEMRXIV_DOI_PREFIX = "10.26434"

# ChemRxiv subject categories (hardcoded from Atypon platform facets)
CHEMRXIV_CATEGORIES: dict[int, str] = {
    502556: "Analytical Chemistry",
    502557: "Biological and Medicinal Chemistry",
    502558: "Catalysis",
    502559: "Chemical Biology",
    502560: "Chemical Engineering and Industrial Chemistry",
    502561: "Earth, Space, and Environmental Chemistry",
    502562: "Education",
    502563: "Inorganic Chemistry",
    502564: "Materials Chemistry",
    502565: "Materials Science",
    502566: "Nanoscience",
    502567: "Organic Chemistry",
    502568: "Organometallic Chemistry",
    502569: "Physical Chemistry",
    502570: "Polymer Chemistry",
    502571: "Supramolecular Chemistry",
    502572: "Theoretical and Computational Chemistry",
    502573: "Other",
}


# =============================================================================
# Shared HTTP helpers
# =============================================================================


@contextlib.asynccontextmanager
async def _http(**kwargs):
    """Shared async HTTP client context manager."""
    defaults = {
        "timeout": TIMEOUT,
        "follow_redirects": True,
        "headers": {"User-Agent": USER_AGENT, "Accept": "application/json"},
    }
    defaults.update(kwargs)
    async with httpx.AsyncClient(**defaults) as client:
        yield client


async def _get(
    url: str,
    params: dict | None = None,
    headers: dict | None = None,
) -> dict | None:
    """HTTP GET returning parsed JSON or None on failure."""
    try:
        async with _http() as client:
            resp = await client.get(url, params=params, headers=headers or {})
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPStatusError as e:
        logger.warning(f"HTTP {e.response.status_code} for GET {url}")
        return None
    except Exception as e:
        logger.warning(f"GET {url} failed: {e}")
        return None


async def _post(
    url: str,
    json_data: dict | None = None,
    params: dict | None = None,
    headers: dict | None = None,
) -> dict | None:
    """HTTP POST returning parsed JSON or None on failure."""
    try:
        async with _http() as client:
            resp = await client.post(
                url, json=json_data, params=params, headers=headers or {}
            )
            resp.raise_for_status()
            return resp.json()
    except httpx.HTTPStatusError as e:
        logger.warning(f"HTTP {e.response.status_code} for POST {url}")
        return None
    except Exception as e:
        logger.warning(f"POST {url} failed: {e}")
        return None


# =============================================================================
# Crossref API
# =============================================================================


async def crossref_search(
    query: str,
    rows: int = 10,
    offset: int = 0,
    sort: str = "relevance",
    filters: dict[str, str] | None = None,
) -> dict | None:
    """Search Crossref works. Returns full API response dict."""
    params: dict[str, Any] = {"query": query, "rows": rows, "offset": offset}
    if sort == "date":
        params["sort"] = "deposited"
        params["order"] = "desc"
    if filters:
        params["filter"] = ",".join(f"{k}:{v}" for k, v in filters.items())
    headers = {}
    if OPENALEX_EMAIL:
        headers["mailto"] = OPENALEX_EMAIL  # Crossref polite pool
    return await _get(f"{CROSSREF_BASE}/works", params=params, headers=headers)


async def crossref_search_chemrxiv(
    query: str = "",
    rows: int = 10,
    offset: int = 0,
    sort: str = "relevance",
    date_from: str | None = None,
    date_to: str | None = None,
) -> dict | None:
    """Search specifically within ChemRxiv preprints via Crossref DOI prefix."""
    filters: dict[str, str] = {"prefix": CHEMRXIV_DOI_PREFIX}
    if date_from:
        filters["from-posted-date"] = date_from
    if date_to:
        filters["until-posted-date"] = date_to
    return await crossref_search(
        query=query, rows=rows, offset=offset, sort=sort, filters=filters
    )


async def crossref_get_work(doi: str) -> dict | None:
    """Get a single work by DOI. Returns the message object."""
    headers = {}
    if OPENALEX_EMAIL:
        headers["mailto"] = OPENALEX_EMAIL
    data = await _get(f"{CROSSREF_BASE}/works/{doi}", headers=headers)
    if data:
        return data.get("message")
    return None


# =============================================================================
# OpenAlex API
# =============================================================================


def _oa_params(**extra) -> dict:
    """Build OpenAlex params with polite pool email."""
    params = dict(extra)
    if OPENALEX_EMAIL:
        params["mailto"] = OPENALEX_EMAIL
    return params


async def openalex_search(
    query: str,
    filters: str | None = None,
    per_page: int = 10,
    page: int = 1,
    sort: str | None = None,
) -> dict | None:
    """Search OpenAlex works. Supports advanced filters and sorting."""
    params = _oa_params(search=query, per_page=per_page, page=page)
    if filters:
        params["filter"] = filters
    if sort:
        params["sort"] = sort
    return await _get(f"{OPENALEX_BASE}/works", params=params)


async def openalex_get_work(identifier: str) -> dict | None:
    """Get work by DOI or OpenAlex ID (W-prefixed)."""
    if identifier.startswith("W") or identifier.startswith("https://"):
        url = f"{OPENALEX_BASE}/works/{identifier}"
    else:
        url = f"{OPENALEX_BASE}/works/doi:{identifier}"
    return await _get(url, params=_oa_params())


async def openalex_get_author(identifier: str) -> dict | None:
    """Get author by OpenAlex ID or search by name."""
    if identifier.startswith("A") or identifier.startswith("https://"):
        return await _get(
            f"{OPENALEX_BASE}/authors/{identifier}", params=_oa_params()
        )
    # Name search — return first result
    data = await _get(
        f"{OPENALEX_BASE}/authors",
        params=_oa_params(search=identifier, per_page=1),
    )
    if data and data.get("results"):
        return data["results"][0]
    return data


async def openalex_search_authors(
    query: str, per_page: int = 10
) -> dict | None:
    """Search authors by name."""
    return await _get(
        f"{OPENALEX_BASE}/authors",
        params=_oa_params(search=query, per_page=per_page),
    )


async def openalex_group_by(
    group_by: str,
    filters: str | None = None,
    search: str | None = None,
    per_page: int = 200,
) -> dict | None:
    """Aggregate works by a field for bibliometric analysis.

    group_by options: publication_year, authorships.author.id,
    primary_location.source.id, topics.id, open_access.oa_status, etc.
    """
    params = _oa_params(group_by=group_by, per_page=per_page)
    if filters:
        params["filter"] = filters
    if search:
        params["search"] = search
    return await _get(f"{OPENALEX_BASE}/works", params=params)


async def openalex_get_topic(identifier: str) -> dict | None:
    """Get topic by OpenAlex ID or search by name."""
    if identifier.startswith("T") or identifier.startswith("https://"):
        return await _get(
            f"{OPENALEX_BASE}/topics/{identifier}", params=_oa_params()
        )
    data = await _get(
        f"{OPENALEX_BASE}/topics",
        params=_oa_params(search=identifier, per_page=1),
    )
    if data and data.get("results"):
        return data["results"][0]
    return data


async def openalex_get_works(
    filters: str,
    per_page: int = 10,
    page: int = 1,
    sort: str | None = None,
) -> dict | None:
    """Get works by filter (no search term). For citation network queries."""
    params = _oa_params(filter=filters, per_page=per_page, page=page)
    if sort:
        params["sort"] = sort
    return await _get(f"{OPENALEX_BASE}/works", params=params)


def openalex_reconstruct_abstract(inverted_index: dict) -> str:
    """Reconstruct abstract text from OpenAlex inverted index."""
    if not inverted_index:
        return ""
    words: dict[int, str] = {}
    for word, positions in inverted_index.items():
        for pos in positions:
            words[pos] = word
    return " ".join(words[i] for i in sorted(words)) if words else ""


# =============================================================================
# Semantic Scholar API
# =============================================================================


def _s2_headers() -> dict:
    """Semantic Scholar headers with optional API key."""
    h: dict[str, str] = {}
    if S2_API_KEY:
        h["x-api-key"] = S2_API_KEY
    return h


async def s2_search(
    query: str,
    limit: int = 10,
    offset: int = 0,
    fields: str = S2_SEARCH_FIELDS,
    fields_of_study: str | None = None,
    year: str | None = None,
    publication_types: str | None = None,
    open_access_pdf: bool | None = None,
) -> dict | None:
    """Search Semantic Scholar papers.

    Supports boolean queries, exact phrases ("..."), and filters.
    year format: "2020" or "2020-2025" or "2020-"
    """
    params: dict[str, Any] = {
        "query": query,
        "fields": fields,
        "limit": limit,
        "offset": offset,
    }
    if fields_of_study:
        params["fieldsOfStudy"] = fields_of_study
    if year:
        params["year"] = year
    if publication_types:
        params["publicationTypes"] = publication_types
    if open_access_pdf:
        params["openAccessPdf"] = ""
    return await _get(
        f"{S2_BASE}/graph/v1/paper/search",
        params=params,
        headers=_s2_headers(),
    )


async def s2_get_paper(
    paper_id: str, fields: str = S2_DETAIL_FIELDS
) -> dict | None:
    """Get paper by Semantic Scholar ID or external ID.

    paper_id formats:
      - S2 paper ID (40-char hex)
      - DOI:10.1234/xxx
      - ARXIV:2101.12345
      - PMID:12345678
      - CorpusId:12345678
    """
    return await _get(
        f"{S2_BASE}/graph/v1/paper/{paper_id}",
        params={"fields": fields},
        headers=_s2_headers(),
    )


async def s2_get_citations(
    paper_id: str,
    limit: int = 50,
    offset: int = 0,
    fields: str = S2_CITATION_FIELDS,
) -> dict | None:
    """Get papers that cite the given paper.

    Includes isInfluential flag and citation intents.
    """
    return await _get(
        f"{S2_BASE}/graph/v1/paper/{paper_id}/citations",
        params={"fields": fields, "limit": limit, "offset": offset},
        headers=_s2_headers(),
    )


async def s2_get_references(
    paper_id: str,
    limit: int = 50,
    offset: int = 0,
    fields: str = S2_CITATION_FIELDS,
) -> dict | None:
    """Get papers referenced by the given paper."""
    return await _get(
        f"{S2_BASE}/graph/v1/paper/{paper_id}/references",
        params={"fields": fields, "limit": limit, "offset": offset},
        headers=_s2_headers(),
    )


async def s2_get_recommendations(
    positive_ids: list[str],
    negative_ids: list[str] | None = None,
    limit: int = 10,
    fields: str = S2_SEARCH_FIELDS,
) -> dict | None:
    """Get paper recommendations based on positive/negative examples.

    IDs can be S2 paper IDs, DOI:xxx, ARXIV:xxx, etc.
    For single-paper similarity, pass one positive ID.
    """
    body: dict[str, Any] = {"positivePaperIds": positive_ids}
    if negative_ids:
        body["negativePaperIds"] = negative_ids
    return await _post(
        f"{S2_BASE}/recommendations/v1/papers/",
        json_data=body,
        params={"fields": fields, "limit": limit},
        headers=_s2_headers(),
    )


async def s2_search_author(
    query: str, limit: int = 5, fields: str = S2_AUTHOR_FIELDS
) -> dict | None:
    """Search for authors by name."""
    return await _get(
        f"{S2_BASE}/graph/v1/author/search",
        params={"query": query, "limit": limit, "fields": fields},
        headers=_s2_headers(),
    )


async def s2_get_author(
    author_id: str, fields: str = S2_AUTHOR_DETAIL_FIELDS
) -> dict | None:
    """Get author details including recent papers."""
    return await _get(
        f"{S2_BASE}/graph/v1/author/{author_id}",
        params={"fields": fields},
        headers=_s2_headers(),
    )


# =============================================================================
# Unpaywall API
# =============================================================================


async def unpaywall_get(doi: str) -> dict | None:
    """Find open access PDF location for a DOI."""
    email = UNPAYWALL_EMAIL or "scholarly-mcp@users.noreply.github.com"
    return await _get(f"{UNPAYWALL_BASE}/{doi}", params={"email": email})


# =============================================================================
# PubChem API
# =============================================================================


async def pubchem_search_by_name(name: str) -> dict | None:
    """Search PubChem compounds by name."""
    return await _get(f"{PUBCHEM_BASE}/compound/name/{name}/JSON")


async def pubchem_search_by_smiles(smiles: str) -> dict | None:
    """Search PubChem compounds by SMILES string."""
    # Use POST for SMILES to handle special characters
    try:
        async with _http() as client:
            resp = await client.post(
                f"{PUBCHEM_BASE}/compound/smiles/JSON",
                data={"smiles": smiles},
            )
            resp.raise_for_status()
            return resp.json()
    except Exception as e:
        logger.warning(f"PubChem SMILES search failed: {e}")
        return None


async def pubchem_search_by_formula(formula: str) -> dict | None:
    """Search PubChem compounds by molecular formula."""
    return await _get(
        f"{PUBCHEM_BASE}/compound/fastformula/{formula}/JSON",
        params={"MaxRecords": 10},
    )


async def pubchem_get_compound(cid: int | str) -> dict | None:
    """Get compound record by PubChem CID."""
    return await _get(f"{PUBCHEM_BASE}/compound/cid/{cid}/JSON")


async def pubchem_get_properties(
    cid: int | str,
    properties: str | None = None,
) -> dict | None:
    """Get computed molecular properties for a compound.

    Default properties cover Lipinski rule-of-5 and common descriptors.
    """
    if properties is None:
        properties = (
            "MolecularFormula,MolecularWeight,CanonicalSMILES,"
            "IsomericSMILES,InChI,InChIKey,IUPACName,XLogP,"
            "ExactMass,MonoisotopicMass,TPSA,Complexity,"
            "HBondDonorCount,HBondAcceptorCount,RotatableBondCount,"
            "HeavyAtomCount,CovalentUnitCount"
        )
    return await _get(
        f"{PUBCHEM_BASE}/compound/cid/{cid}/property/{properties}/JSON"
    )


async def pubchem_get_synonyms(cid: int | str) -> dict | None:
    """Get all known names/synonyms for a compound."""
    return await _get(f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON")


# =============================================================================
# Common Chemistry (CAS) API — free, no auth, ~500k compounds
# =============================================================================


async def cas_search(query: str, size: int = 10) -> dict | None:
    """Search Common Chemistry by name, CAS number, InChI, InChIKey, or SMILES."""
    return await _get(f"{CAS_BASE}/search", params={"q": query, "size": size})


async def cas_detail(cas_rn: str) -> dict | None:
    """Get full compound details by CAS Registry Number.

    Returns: name, CAS RN, molecular formula, molecular mass, InChI,
    InChIKey, SMILES, canonical SMILES, and experimental properties.
    """
    return await _get(f"{CAS_BASE}/detail", params={"cas_rn": cas_rn})


# =============================================================================
# Web of Science Starter API (optional — requires WOS_API_KEY)
# =============================================================================


def wos_available() -> bool:
    """Check if Web of Science API credentials are configured."""
    return bool(WOS_API_KEY)


async def wos_search(
    query: str,
    limit: int = 10,
    first_record: int = 1,
    sort_field: str = "RS",
    database_id: str = "WOS",
) -> dict | None:
    """Search Web of Science.

    Requires WOS_API_KEY environment variable.

    query: Web of Science advanced search syntax, e.g.:
      - TS=(catalysis AND asymmetric)  — topic search
      - AU=(Smith)                     — author search
      - SO=(Nature)                    — source/journal search
      - DO=(10.1234/xxx)               — DOI search

    sort_field: RS (relevance), PY (year), TC (times cited), LD (load date)
    database_id: WOS, BCI, CCC, DCI, DIIDW, KJD, MEDLINE, RSCI, SCIELO
    """
    if not WOS_API_KEY:
        return None
    return await _get(
        WOS_BASE,
        params={
            "databaseId": database_id,
            "usrQuery": query,
            "count": limit,
            "firstRecord": first_record,
            "sortField": sort_field,
        },
        headers={"X-ApiKey": WOS_API_KEY},
    )

# =============================================================================
# NIST Chemistry WebBook (scraping — no official API)
# =============================================================================


def _strip_html(html: str) -> str:
    """Remove HTML tags and decode common entities."""
    text = re.sub(r"<[^>]+>", "", html)
    for ent, char in [("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"),
                       ("&plusmn;", "±"), ("&deg;", "°"), ("&nbsp;", " "),
                       ("&#176;", "°")]:
        text = text.replace(ent, char)
    return text.strip()


async def nist_fetch(params: dict[str, str]) -> str | None:
    """Fetch a NIST WebBook page. Returns raw HTML."""
    params.setdefault("Units", "SI")
    try:
        async with _http() as client:
            resp = await client.get(NIST_BASE, params=params)
            if resp.status_code == 200:
                return resp.text
    except Exception as e:
        logger.warning(f"NIST fetch failed: {e}")
    return None


async def nist_search(query: str, search_type: str = "name") -> str | None:
    """Search NIST WebBook by name, CAS, formula, or InChI."""
    params: dict[str, str] = {}
    if search_type == "cas":
        cas_clean = "C" + query.replace("-", "")
        params["ID"] = cas_clean
        params["Mask"] = "FFF"
    elif search_type == "formula":
        params["Formula"] = query
        params["NoIon"] = "on"
    elif search_type == "inchi":
        params["InChI"] = query
    else:
        params["Name"] = query
    return await nist_fetch(params)


def nist_is_compound_page(html: str) -> bool:
    """Check if HTML is a single compound page (vs. search results)."""
    return bool(re.search(r'<h1[^>]*id="Top"', html))


def nist_parse_search_results(html: str) -> list[dict]:
    """Parse NIST search results page into list of matches."""
    results = []
    for m in re.finditer(
        r'<a\s+href="(/cgi/cbook\.cgi\?ID=(C\d+)[^"]*)"[^>]*>(.*?)</a>',
        html,
    ):
        results.append({
            "name": _strip_html(m.group(3)),
            "nist_id": m.group(2),
            "url": f"https://webbook.nist.gov{m.group(1)}",
        })
    return results


def nist_parse_compound(html: str) -> dict:
    """Parse NIST compound page into structured data."""
    info: dict[str, Any] = {}

    # Name
    m = re.search(r'<h1[^>]*>(.*?)</h1>', html, re.S)
    if m:
        info["name"] = _strip_html(m.group(1))

    # Key-value pairs from <li><strong>Key:</strong> Value</li>
    kv_map = {
        "Formula": "formula",
        "Molecular weight": "molecular_weight",
        "CAS Registry Number": "cas_rn",
        "IUPAC Standard InChI": "inchi",
        "IUPAC Standard InChIKey": "inchi_key",
        "Chemical structure": None,  # skip image
    }
    for m in re.finditer(
        r"<li>\s*<strong>(.*?):</strong>\s*(.*?)</li>", html, re.S
    ):
        key = _strip_html(m.group(1))
        val = _strip_html(m.group(2))
        mapped = kv_map.get(key, key.lower().replace(" ", "_"))
        if mapped and val:
            info[mapped] = val

    # Other names
    m = re.search(r"<strong>Other names:</strong>\s*(.*?)</li>", html, re.S)
    if m:
        raw = _strip_html(m.group(1))
        info["other_names"] = [n.strip() for n in raw.split(";") if n.strip()]

    # NIST ID from page URLs
    m = re.search(r"ID=(C\d+)", html)
    if m:
        info["nist_id"] = m.group(1)

    # Detect available data sections
    avail: list[str] = []
    section_names = [
        ("Thermochemistry", "thermo"), ("Phase change", "phase_change"),
        ("Reaction thermochemistry", "reaction_thermo"),
        ("Henry", "henrys_law"), ("Gas phase ion", "ion_energetics"),
        ("IR Spec", "ir_spectrum"), ("Mass Spec", "mass_spectrum"),
        ("UV/Vis", "uv_vis_spectrum"), ("Vibrational", "vibrational"),
        ("Electronic", "electronic"), ("Constants of diatomic", "diatomic"),
    ]
    for label, key in section_names:
        if label.lower() in html.lower():
            avail.append(key)
    info["available_data"] = avail

    # --- Inline thermochemistry data ---
    # Gas phase ΔfH°
    m = re.search(
        r"f</sub>H.*?gas.*?(-?[\d.]+)\s*(?:±|&plusmn;)\s*([\d.]+)\s*kJ/mol",
        html, re.S,
    )
    if m:
        info["delta_fH_gas_kJ_mol"] = f"{m.group(1)} ± {m.group(2)}"

    # Standard entropy S°
    m = re.search(r"S°.*?gas.*?([\d.]+)\s*(?:±|&plusmn;)?\s*[\d.]*\s*J/mol", html, re.S)
    if m:
        info["S_gas_J_mol_K"] = m.group(1)

    # Cp gas
    m = re.search(r"C\s*p.*?gas.*?([\d.]+)\s*(?:±|&plusmn;)?\s*[\d.]*\s*J/mol", html, re.S)
    if m:
        info["Cp_gas_J_mol_K"] = m.group(1)

    # Phase change: boiling point
    for pat in [
        r"T<sub>boil</sub>\s*=?\s*([\d.]+)\s*(?:±\s*[\d.]+\s*)?K",
        r"boil.*?([\d.]+)\s*K",
    ]:
        m = re.search(pat, html, re.S)
        if m:
            info["boiling_point_K"] = m.group(1)
            break

    # Phase change: melting point
    for pat in [
        r"T<sub>fus</sub>\s*=?\s*([\d.]+)\s*(?:±\s*[\d.]+\s*)?K",
        r"fus.*?([\d.]+)\s*K",
    ]:
        m = re.search(pat, html, re.S)
        if m:
            info["melting_point_K"] = m.group(1)
            break

    # ΔvapH (enthalpy of vaporization)
    m = re.search(r"vap</sub>H.*?([\d.]+)\s*(?:±|&plusmn;)\s*([\d.]+)\s*kJ/mol", html, re.S)
    if m:
        info["delta_vapH_kJ_mol"] = f"{m.group(1)} ± {m.group(2)}"

    # ΔfusH (enthalpy of fusion)
    m = re.search(r"fus</sub>H.*?([\d.]+)\s*(?:±|&plusmn;)\s*([\d.]+)\s*kJ/mol", html, re.S)
    if m:
        info["delta_fusH_kJ_mol"] = f"{m.group(1)} ± {m.group(2)}"

    return info


# =============================================================================
# Materials Project API (optional — requires MP_API_KEY)
# =============================================================================


def mp_available() -> bool:
    """Check if Materials Project API key is configured."""
    return bool(MP_API_KEY)


async def mp_search(
    formula: str | None = None,
    elements: list[str] | None = None,
    band_gap_min: float | None = None,
    band_gap_max: float | None = None,
    limit: int = 10,
) -> dict | None:
    """Search Materials Project for inorganic materials."""
    if not MP_API_KEY:
        return None
    params: dict[str, Any] = {
        "_limit": limit,
        "_fields": (
            "material_id,formula_pretty,structure,symmetry,"
            "band_gap,formation_energy_per_atom,energy_above_hull,"
            "is_stable,theoretical,nsites"
        ),
    }
    if formula:
        params["formula"] = formula
    if elements:
        params["elements"] = ",".join(elements)
    if band_gap_min is not None:
        params["band_gap_min"] = band_gap_min
    if band_gap_max is not None:
        params["band_gap_max"] = band_gap_max
    return await _get(
        f"{MP_BASE}/materials/summary/",
        params=params,
        headers={"X-API-KEY": MP_API_KEY},
    )


async def mp_get_material(material_id: str) -> dict | None:
    """Get full material details by Materials Project ID (e.g., 'mp-149')."""
    if not MP_API_KEY:
        return None
    return await _get(
        f"{MP_BASE}/materials/summary/{material_id}",
        params={
            "_fields": (
                "material_id,formula_pretty,structure,symmetry,"
                "band_gap,formation_energy_per_atom,energy_above_hull,"
                "is_stable,theoretical,nsites,volume,density,"
                "efermi,total_magnetization,ordering,is_metal,"
                "database_IDs,deprecated,uncorrected_energy_per_atom"
            ),
        },
        headers={"X-API-KEY": MP_API_KEY},
    )


# =============================================================================
# RXN4Chemistry — IBM AI reaction prediction (optional — requires RXN_API_KEY)
# =============================================================================


def rxn_available() -> bool:
    """Check if IBM RXN API key is configured."""
    return bool(RXN_API_KEY)


_rxn_project_id: str | None = None


async def _rxn_headers() -> dict[str, str]:
    """Get RXN auth headers."""
    return {
        "Authorization": f"apikey {RXN_API_KEY}",
        "Content-Type": "application/json",
    }


async def _rxn_ensure_project() -> str | None:
    """Create or return cached RXN project ID."""
    global _rxn_project_id
    if _rxn_project_id:
        return _rxn_project_id
    if not RXN_API_KEY:
        return None
    data = await _post(
        f"{RXN_BASE}/projects",
        json_data={"name": "scholarly-mcp"},
        headers=await _rxn_headers(),
    )
    if data and data.get("payload"):
        _rxn_project_id = data["payload"].get("id")
    return _rxn_project_id


async def _rxn_poll(url: str, max_wait: int = 60, interval: int = 3) -> dict | None:
    """Poll an RXN endpoint until result is ready."""
    headers = await _rxn_headers()
    for _ in range(max_wait // interval):
        data = await _get(url, headers=headers)
        if data:
            payload = data.get("payload", {})
            status = payload.get("status", data.get("status", ""))
            if status in ("SUCCESS", "success"):
                return data
            if status in ("FAILED", "failed", "ERROR", "error"):
                return data
        await asyncio.sleep(interval)
    return None


async def rxn_predict_reaction(rxn_smiles: str) -> dict | None:
    """Predict reaction product(s) from reactants.

    rxn_smiles: reactants.reagents>>products (e.g., 'CC.OC>>')
    """
    if not RXN_API_KEY:
        return None
    project_id = await _rxn_ensure_project()
    if not project_id:
        return None
    data = await _post(
        f"{RXN_BASE}/predictions",
        json_data={
            "projectId": project_id,
            "name": "mcp-prediction",
            "inputs": [{"rxnSmiles": rxn_smiles}],
        },
        headers=await _rxn_headers(),
    )
    if not data or not data.get("payload"):
        return data
    pred_id = data["payload"].get("id")
    if not pred_id:
        return data
    return await _rxn_poll(
        f"{RXN_BASE}/predictions/{pred_id}",
    )


async def rxn_retrosynthesis(product_smiles: str, max_steps: int = 3) -> dict | None:
    """Plan retrosynthetic route to a target molecule.

    product_smiles: target product SMILES
    max_steps: maximum retrosynthetic steps (1-10)
    """
    if not RXN_API_KEY:
        return None
    project_id = await _rxn_ensure_project()
    if not project_id:
        return None
    data = await _post(
        f"{RXN_BASE}/retrosynthesis",
        json_data={
            "projectId": project_id,
            "fap": 0.6,
            "maxSteps": max_steps,
            "nBeams": 10,
            "pruneThreshold": 0.2,
            "isAutomatic": True,
            "content": {"smiles": product_smiles},
        },
        headers=await _rxn_headers(),
    )
    if not data or not data.get("payload"):
        return data
    pred_id = data["payload"].get("id")
    if not pred_id:
        return data
    # Retrosynthesis takes longer — poll with longer timeout
    return await _rxn_poll(
        f"{RXN_BASE}/retrosynthesis/{pred_id}?projectId={project_id}",
        max_wait=120,
        interval=5,
    )


# =============================================================================
# UniChem — Universal chemical identifier cross-reference (no auth)
# =============================================================================


async def unichem_lookup(inchikey: str) -> dict | None:
    """Cross-reference a compound across 40+ databases by InChIKey.

    Returns source IDs from ChEMBL, PubChem, DrugBank, ZINC, etc.
    """
    return await _post(
        f"{UNICHEM_BASE}/compounds",
        json_data={"type": "inchikey", "compound": inchikey},
    )


async def unichem_sources() -> dict | None:
    """List all available UniChem data sources."""
    return await _get(f"{UNICHEM_BASE}/sources")


# =============================================================================
# Crystallography Open Database (COD) — open crystal structures (no auth)
# =============================================================================


async def cod_search(
    formula: str | None = None,
    elements: list[str] | None = None,
    text: str | None = None,
    limit: int = 20,
) -> list | None:
    """Search COD for crystal structures.

    formula: Hill notation (e.g., 'C6 H6', 'Fe2 O3')
    elements: required elements (e.g., ['Fe', 'O'])
    text: free text search in compound names
    """
    params: dict[str, str] = {"format": "json"}
    if formula:
        params["formula"] = formula
    if elements:
        for i, el in enumerate(elements[:8], 1):
            params[f"el{i}"] = el
    if text:
        params["text"] = text
    try:
        async with _http() as client:
            resp = await client.get(f"{COD_BASE}/result", params=params)
            if resp.status_code == 200:
                data = resp.json()
                if isinstance(data, list):
                    return data[:limit]
                return data
    except Exception as e:
        logger.warning(f"COD search failed: {e}")
    return None


async def cod_get_cif(cod_id: int | str) -> str | None:
    """Download CIF file for a COD entry."""
    try:
        async with _http() as client:
            resp = await client.get(f"{COD_BASE}/{cod_id}.cif")
            if resp.status_code == 200:
                return resp.text
    except Exception as e:
        logger.warning(f"COD CIF download failed: {e}")
    return None


# =============================================================================
# EPA CompTox Dashboard (optional — requires COMPTOX_API_KEY)
# =============================================================================


def comptox_available() -> bool:
    """Check if EPA CompTox API key is configured."""
    return bool(COMPTOX_API_KEY)


async def comptox_search(query: str) -> dict | None:
    """Search CompTox by chemical name, CAS, or DTXSID."""
    if not COMPTOX_API_KEY:
        return None
    # Try name search
    return await _get(
        f"{COMPTOX_BASE}/chemical/search/by-name/{query}",
        headers={"x-api-key": COMPTOX_API_KEY, "Accept": "application/json"},
    )


async def comptox_get_details(dtxsid: str) -> dict | None:
    """Get full chemical details by DTXSID identifier."""
    if not COMPTOX_API_KEY:
        return None
    return await _get(
        f"{COMPTOX_BASE}/chemical/detail/search/by-dtxsid/{dtxsid}",
        headers={"x-api-key": COMPTOX_API_KEY, "Accept": "application/json"},
    )


async def comptox_get_properties(dtxsid: str) -> dict | None:
    """Get physicochemical and fate properties for a chemical."""
    if not COMPTOX_API_KEY:
        return None
    return await _get(
        f"{COMPTOX_BASE}/chemical/property/search/by-dtxsid/{dtxsid}",
        headers={"x-api-key": COMPTOX_API_KEY, "Accept": "application/json"},
    )


async def comptox_get_hazard(dtxsid: str) -> dict | None:
    """Get hazard data for a chemical."""
    if not COMPTOX_API_KEY:
        return None
    return await _get(
        f"{COMPTOX_BASE}/hazard/search/by-dtxsid/{dtxsid}",
        headers={"x-api-key": COMPTOX_API_KEY, "Accept": "application/json"},
    )


# =============================================================================
# MassBank EU — Mass spectrometry reference spectra (no auth)
# =============================================================================


async def massbank_search(
    compound_name: str | None = None,
    formula: str | None = None,
    inchikey: str | None = None,
    exact_mass_min: float | None = None,
    exact_mass_max: float | None = None,
    instrument_type: str | None = None,
    limit: int = 20,
) -> list | None:
    """Search MassBank for reference mass spectra."""
    params: dict[str, Any] = {"limit": limit}
    if compound_name:
        params["compound_name"] = compound_name
    if formula:
        params["formula"] = formula
    if inchikey:
        params["inchi_key"] = inchikey
    if exact_mass_min is not None:
        params["exact_mass_from"] = exact_mass_min
    if exact_mass_max is not None:
        params["exact_mass_to"] = exact_mass_max
    if instrument_type:
        params["instrument_type"] = instrument_type
    try:
        async with _http() as client:
            resp = await client.get(f"{MASSBANK_BASE}/records", params=params)
            if resp.status_code == 200:
                data = resp.json()
                return data if isinstance(data, list) else data.get("data", [])
    except Exception as e:
        logger.warning(f"MassBank search failed: {e}")
    return None


async def massbank_get_record(accession: str) -> dict | None:
    """Get a specific MassBank spectrum record."""
    return await _get(f"{MASSBANK_BASE}/records/{accession}")


# =============================================================================
# BindingDB — Protein-ligand binding affinities (no auth)
# =============================================================================


async def bindingdb_by_target(
    uniprot_id: str,
    cutoff_nm: int = 10000,
) -> list[dict] | None:
    """Get ligands for a protein target by UniProt ID.

    cutoff_nm: binding affinity cutoff in nM (default 10 µM)
    Returns list of dicts with compound SMILES, Ki, IC50, Kd, EC50.
    """
    try:
        async with _http() as client:
            resp = await client.get(
                f"{BINDINGDB_BASE}/getLigandsByUniprots",
                params={
                    "uniprot": uniprot_id,
                    "cutoff": cutoff_nm,
                    "response": "application/json",
                },
                timeout=60,
            )
            if resp.status_code == 200:
                # BindingDB may return JSON or TSV depending on version
                try:
                    return resp.json()
                except Exception:
                    # Parse TSV fallback
                    return _parse_bindingdb_tsv(resp.text)
    except Exception as e:
        logger.warning(f"BindingDB target search failed: {e}")
    return None


async def bindingdb_by_smiles(
    smiles: str,
    cutoff: float = 0.8,
) -> list[dict] | None:
    """Find similar compounds in BindingDB by SMILES.

    cutoff: Tanimoto similarity threshold (0-1, default 0.8)
    """
    try:
        async with _http() as client:
            resp = await client.get(
                f"{BINDINGDB_BASE}/getTargetByCompound",
                params={
                    "smiles": smiles,
                    "cutoff": cutoff,
                    "response": "application/json",
                },
                timeout=60,
            )
            if resp.status_code == 200:
                try:
                    return resp.json()
                except Exception:
                    return _parse_bindingdb_tsv(resp.text)
    except Exception as e:
        logger.warning(f"BindingDB SMILES search failed: {e}")
    return None


def _parse_bindingdb_tsv(text: str) -> list[dict]:
    """Parse BindingDB tab-separated response into list of dicts."""
    lines = text.strip().split("\n")
    if len(lines) < 2:
        return []
    headers = lines[0].split("\t")
    results = []
    for line in lines[1:]:
        vals = line.split("\t")
        row = {}
        for i, h in enumerate(headers):
            if i < len(vals):
                row[h.strip()] = vals[i].strip()
        results.append(row)
    return results


# =============================================================================
# Crossref BibTeX (content negotiation — no extra API)
# =============================================================================


async def crossref_get_bibtex(doi: str) -> str | None:
    """Get BibTeX entry for a DOI via Crossref content negotiation."""
    doi = doi.strip().removeprefix("https://doi.org/").removeprefix("http://doi.org/")
    try:
        async with _http() as client:
            resp = await client.get(
                f"https://doi.org/{doi}",
                headers={"Accept": "application/x-bibtex"},
                follow_redirects=True,
            )
            if resp.status_code == 200 and "@" in resp.text:
                return resp.text.strip()
    except Exception as e:
        logger.warning(f"BibTeX fetch failed for {doi}: {e}")
    return None


async def crossref_get_bibtex_batch(dois: list[str]) -> list[tuple[str, str | None]]:
    """Get BibTeX entries for multiple DOIs. Returns list of (doi, bibtex)."""
    results = []
    for doi in dois:
        bib = await crossref_get_bibtex(doi)
        results.append((doi, bib))
    return results


# =============================================================================
# RCSB PDB — Protein Data Bank (no auth)
# =============================================================================


async def pdb_search(
    query: str,
    search_type: str = "full_text",
    limit: int = 10,
) -> dict | None:
    """Search RCSB PDB for protein/nucleic acid structures.

    search_type: 'full_text', 'structure_title', 'structure_author'
    """
    service_map = {
        "full_text": "full_text",
        "structure_title": "text",
        "structure_author": "text",
    }
    service = service_map.get(search_type, "full_text")

    json_body: dict[str, Any] = {
        "query": {
            "type": "terminal",
            "service": service,
            "parameters": {"value": query},
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "paginate": {"start": 0, "rows": limit},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }

    # For author/title, use the text service with specific attribute
    if search_type == "structure_title":
        json_body["query"]["parameters"] = {
            "attribute": "struct.title",
            "operator": "contains_phrase",
            "value": query,
        }
    elif search_type == "structure_author":
        json_body["query"]["parameters"] = {
            "attribute": "rcsb_primary_citation.rcsb_authors",
            "operator": "contains_phrase",
            "value": query,
        }

    try:
        async with _http() as client:
            resp = await client.post(
                PDB_SEARCH_BASE,
                json=json_body,
                timeout=30,
            )
            if resp.status_code == 200:
                return resp.json()
    except Exception as e:
        logger.warning(f"PDB search failed: {e}")
    return None


async def pdb_get_entry(pdb_id: str) -> dict | None:
    """Get full entry details from RCSB PDB."""
    pdb_id = pdb_id.strip().upper()
    return await _get(f"{PDB_DATA_BASE}/entry/{pdb_id}")


async def pdb_get_entity(pdb_id: str, entity_id: int = 1) -> dict | None:
    """Get polymer entity details (protein/nucleic acid chain)."""
    pdb_id = pdb_id.strip().upper()
    return await _get(f"{PDB_DATA_BASE}/polymer_entity/{pdb_id}/{entity_id}")


async def pdb_get_ligands(pdb_id: str) -> list[dict]:
    """Get all non-polymer (ligand) entities in a PDB structure."""
    pdb_id = pdb_id.strip().upper()
    ligands = []
    # PDB structures can have multiple non-polymer entities
    for entity_id in range(1, 20):  # usually < 10
        data = await _get(
            f"{PDB_DATA_BASE}/nonpolymer_entity/{pdb_id}/{entity_id}"
        )
        if data:
            ligands.append(data)
        else:
            break
    return ligands


# =============================================================================
# PubChem GHS Hazard Data (extends existing PubChem client)
# =============================================================================


async def pubchem_get_ghs(cid: int) -> dict | None:
    """Get GHS Classification data (hazard pictograms, H/P statements) for a compound.

    Uses PubChem PUG-View API to get the GHS section.
    """
    try:
        async with _http() as client:
            resp = await client.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON",
                params={"heading": "GHS Classification"},
                timeout=20,
            )
            if resp.status_code == 200:
                return resp.json()
    except Exception as e:
        logger.warning(f"PubChem GHS fetch failed for CID {cid}: {e}")
    return None


def parse_ghs_data(pug_view_data: dict) -> dict:
    """Parse PubChem PUG-View GHS response into structured hazard data."""
    result: dict[str, Any] = {
        "pictograms": [],
        "signal_word": "",
        "hazard_statements": [],
        "precautionary_statements": [],
    }

    if not pug_view_data:
        return result

    # Navigate the nested PUG-View structure
    record = pug_view_data.get("Record", {})
    sections = record.get("Section", [])

    for section in sections:
        for subsec in section.get("Section", []):
            heading = subsec.get("TOCHeading", "")

            for info in subsec.get("Information", []):
                val = info.get("Value", {})

                if "Pictogram" in heading or "Pictogram" in info.get("Name", ""):
                    # Extract pictogram names
                    for sv in val.get("StringWithMarkup", []):
                        text = sv.get("String", "")
                        if text:
                            result["pictograms"].append(text)
                        # Also check for markup references
                        for mu in sv.get("Markup", []):
                            extra = mu.get("Extra", "")
                            if extra:
                                result["pictograms"].append(extra)

                elif "Signal" in heading or "Signal" in info.get("Name", ""):
                    for sv in val.get("StringWithMarkup", []):
                        text = sv.get("String", "")
                        if text and text.lower() in ("danger", "warning"):
                            result["signal_word"] = text

                elif "Hazard Statement" in heading or "H Statement" in info.get("Name", ""):
                    for sv in val.get("StringWithMarkup", []):
                        text = sv.get("String", "")
                        if text:
                            result["hazard_statements"].append(text)

                elif "Precautionary" in heading or "P Statement" in info.get("Name", ""):
                    for sv in val.get("StringWithMarkup", []):
                        text = sv.get("String", "")
                        if text:
                            result["precautionary_statements"].append(text)

    # Deduplicate
    result["pictograms"] = list(dict.fromkeys(result["pictograms"]))
    result["hazard_statements"] = list(dict.fromkeys(result["hazard_statements"]))
    result["precautionary_statements"] = list(dict.fromkeys(result["precautionary_statements"]))

    return result


# =============================================================================
# GNPS NPClassifier — Natural product classification (no auth)
# =============================================================================


async def gnps_classify_compound(smiles: str) -> dict | None:
    """Classify a compound into natural product classes using GNPS NPClassifier.

    Returns pathway, superclass, class, and isglycoside prediction.
    """
    try:
        async with _http() as client:
            resp = await client.get(
                f"{NPCLASSIFIER_BASE}/classify",
                params={"smiles": smiles},
                timeout=30,
            )
            if resp.status_code == 200:
                return resp.json()
    except Exception as e:
        logger.warning(f"NPClassifier failed: {e}")
    return None


# =============================================================================
# OpenAlex Sources — Journal metrics (extends existing OpenAlex client)
# =============================================================================


async def openalex_get_source(source_id: str) -> dict | None:
    """Get journal/source details from OpenAlex.

    source_id: OpenAlex source ID (e.g., 'S137773608') or ISSN
    """
    headers = _oa_headers()
    # Try direct ID lookup
    if source_id.startswith("S") or source_id.startswith("https://"):
        return await _get(
            f"https://api.openalex.org/sources/{source_id}",
            headers=headers,
        )
    # Try ISSN lookup
    return await _get(
        f"https://api.openalex.org/sources/issn:{source_id}",
        headers=headers,
    )


async def openalex_search_sources(
    query: str,
    limit: int = 10,
) -> dict | None:
    """Search for journals/sources by name in OpenAlex."""
    return await _get(
        "https://api.openalex.org/sources",
        params={"search": query, "per_page": limit},
        headers=_oa_headers(),
    )
