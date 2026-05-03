"""
NCBI MCP Server — wraps NCBI E-utilities via Biopython Entrez.
Exposes tools for searching databases, fetching records, taxonomy lookup,
genome assembly info, and PubMed abstract retrieval.

Rate limits: 3 req/sec without API key, 10/sec with one.
Set NCBI_API_KEY and NCBI_EMAIL in environment before starting.
"""

import os
import time
import json
import textwrap
from io import StringIO

from mcp.server.fastmcp import FastMCP
from Bio import Entrez, SeqIO
from Bio.Entrez import efetch, esearch, esummary, einfo, elink

# --- Configuration -----------------------------------------------------------

Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY", None)

# 0.11s gap → ~9 req/sec with key; 0.34s → ~3/sec without
_RATE_DELAY = 0.11 if Entrez.api_key else 0.34

mcp = FastMCP("ncbi-server")


def _rate_limit():
    time.sleep(_RATE_DELAY)


# --- Tools -------------------------------------------------------------------

@mcp.tool()
def ncbi_search(db: str, term: str, retmax: int = 20) -> str:
    """
    Search any NCBI database (pubmed, nucleotide, protein, genome, assembly,
    taxonomy, gene, sra, etc.) and return a list of matching IDs.

    Args:
        db:     NCBI database name (e.g. 'pubmed', 'nucleotide', 'assembly').
        term:   Entrez query string (supports field tags like [Organism], [Title]).
        retmax: Maximum number of IDs to return (default 20, max 10000).

    Returns:
        JSON with 'count' (total hits) and 'ids' (list of matching IDs).
    """
    _rate_limit()
    handle = esearch(db=db, term=term, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return json.dumps({
        "count": int(record["Count"]),
        "ids": record["IdList"],
        "query_translation": record.get("QueryTranslation", ""),
    }, indent=2)


@mcp.tool()
def ncbi_fetch(db: str, ids: str, rettype: str = "gb", retmode: str = "text") -> str:
    """
    Fetch full records from NCBI by ID or comma-separated list of IDs.

    Args:
        db:      NCBI database (e.g. 'nucleotide', 'protein', 'pubmed').
        ids:     Comma-separated accession numbers or UIDs.
        rettype: Return type — 'gb' (GenBank), 'fasta', 'abstract', 'medline', etc.
        retmode: Return mode — 'text' or 'xml'.

    Returns:
        Raw record text (GenBank, FASTA, abstract, etc.).
    """
    _rate_limit()
    handle = efetch(db=db, id=ids, rettype=rettype, retmode=retmode)
    data = handle.read()
    handle.close()
    return data if isinstance(data, str) else data.decode("utf-8")


@mcp.tool()
def ncbi_summary(db: str, ids: str) -> str:
    """
    Retrieve document summaries (DocSums) for one or more NCBI records.
    Lighter than ncbi_fetch — useful for metadata without downloading full records.

    Args:
        db:  NCBI database name.
        ids: Comma-separated UIDs or accession numbers.

    Returns:
        JSON-formatted list of document summary dicts.
    """
    _rate_limit()
    handle = esummary(db=db, id=ids, retmode="json")
    record = Entrez.read(handle)
    handle.close()
    return json.dumps(record, indent=2, default=str)


@mcp.tool()
def pubmed_fetch_abstract(pmid: str) -> str:
    """
    Fetch the full abstract and metadata for a PubMed article by PMID.

    Args:
        pmid: PubMed ID (numeric string).

    Returns:
        Structured text with Title, Authors, Journal, Year, and Abstract.
    """
    _rate_limit()
    handle = efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
    text = handle.read()
    handle.close()
    return text


@mcp.tool()
def ncbi_taxonomy(name: str) -> str:
    """
    Resolve an organism name to its NCBI Taxonomy ID and full lineage.
    Always call this before fetching genomes to confirm the correct taxon.

    Args:
        name: Scientific name or common name (e.g. 'Klebsiella pneumoniae').

    Returns:
        JSON with TaxId, ScientificName, Rank, and Lineage.
    """
    _rate_limit()
    search_handle = esearch(db="taxonomy", term=name)
    search_record = Entrez.read(search_handle)
    search_handle.close()

    if not search_record["IdList"]:
        return json.dumps({"error": f"No taxonomy record found for '{name}'"})

    taxid = search_record["IdList"][0]
    _rate_limit()
    fetch_handle = efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(fetch_handle)
    fetch_handle.close()

    rec = records[0]
    return json.dumps({
        "TaxId": rec["TaxId"],
        "ScientificName": rec["ScientificName"],
        "Rank": rec["Rank"],
        "Lineage": rec["Lineage"],
        "GeneticCode": rec.get("GeneticCode", {}).get("GCName", ""),
    }, indent=2)


@mcp.tool()
def ncbi_genome_assembly(organism: str, retmax: int = 5) -> str:
    """
    Search for genome assemblies for a given organism in NCBI Assembly database.
    Returns accession numbers, assembly level, and submission dates.

    Args:
        organism: Organism name (e.g. 'Escherichia coli O157').
        retmax:   Max number of assemblies to return (default 5).

    Returns:
        JSON list of assembly summaries with accession and quality metrics.
    """
    _rate_limit()
    term = f"{organism}[Organism] AND latest[filter]"
    search_handle = esearch(db="assembly", term=term, retmax=retmax)
    search_record = Entrez.read(search_handle)
    search_handle.close()

    if not search_record["IdList"]:
        return json.dumps({"error": f"No assemblies found for '{organism}'"})

    _rate_limit()
    summary_handle = esummary(db="assembly", id=",".join(search_record["IdList"]), retmode="json")
    summaries = Entrez.read(summary_handle)
    summary_handle.close()

    results = []
    for uid, doc in summaries["DocumentSummarySet"]["DocumentSummary"].items():
        results.append({
            "uid": uid,
            "AssemblyAccession": doc.get("AssemblyAccession"),
            "AssemblyName": doc.get("AssemblyName"),
            "AssemblyLevel": doc.get("AssemblyStatus"),
            "Organism": doc.get("Organism"),
            "SubmissionDate": doc.get("SubmissionDate"),
            "ContigN50": doc.get("ContigN50"),
            "ScaffoldN50": doc.get("ScaffoldN50"),
            "Coverage": doc.get("Coverage"),
        })
    return json.dumps(results, indent=2, default=str)


@mcp.tool()
def ncbi_pubmed_search(query: str, retmax: int = 10) -> str:
    """
    Search PubMed and return titles, PMIDs, and journal info for matching articles.
    Use field tags for precision: [Title/Abstract], [Author], [MeSH Terms], [PDAT].

    Args:
        query:  PubMed query string.
        retmax: Number of results to return (default 10).

    Returns:
        JSON list of articles with PMID, title, authors, journal, and year.
    """
    _rate_limit()
    search_handle = esearch(db="pubmed", term=query, retmax=retmax, sort="relevance")
    search_record = Entrez.read(search_handle)
    search_handle.close()

    if not search_record["IdList"]:
        return json.dumps({"count": 0, "articles": []})

    _rate_limit()
    summary_handle = esummary(db="pubmed", id=",".join(search_record["IdList"]), retmode="json")
    summaries = Entrez.read(summary_handle)
    summary_handle.close()

    articles = []
    for uid in search_record["IdList"]:
        doc = summaries["result"].get(uid, {})
        articles.append({
            "pmid": uid,
            "title": doc.get("title", ""),
            "authors": [a.get("name", "") for a in doc.get("authors", [])],
            "journal": doc.get("fulljournalname", ""),
            "year": doc.get("pubdate", "")[:4],
            "doi": next((a["value"] for a in doc.get("articleids", []) if a["idtype"] == "doi"), None),
        })
    return json.dumps({"count": int(search_record["Count"]), "articles": articles}, indent=2)


@mcp.tool()
def ncbi_elink(db_from: str, db_to: str, ids: str) -> str:
    """
    Find linked records across NCBI databases (e.g. gene → protein, pubmed → pmc).

    Args:
        db_from: Source database (e.g. 'gene').
        db_to:   Target database (e.g. 'protein').
        ids:     Comma-separated UIDs in the source database.

    Returns:
        JSON with linked IDs in the target database.
    """
    _rate_limit()
    handle = elink(dbfrom=db_from, db=db_to, id=ids)
    record = Entrez.read(handle)
    handle.close()

    linked_ids = []
    for linkset in record:
        for link_db in linkset.get("LinkSetDb", []):
            if link_db["DbTo"] == db_to:
                linked_ids.extend([lnk["Id"] for lnk in link_db["Link"]])

    return json.dumps({"db_from": db_from, "db_to": db_to, "linked_ids": linked_ids}, indent=2)


# --- Entry point -------------------------------------------------------------

if __name__ == "__main__":
    if not Entrez.email:
        raise EnvironmentError(
            "NCBI_EMAIL environment variable is required. "
            "Set it before starting the server: export NCBI_EMAIL=you@example.com"
        )
    mcp.run(transport="stdio")
