from Bio import Entrez, SeqIO
from config.settings import configure_entrez as _central_configure_entrez
from io import StringIO
import json
import re


def configure_entrez(email: str = "", api_key: str = ""):
    """Backward-compatible wrapper around the central NCBI configuration."""
    _central_configure_entrez(email=email or None, api_key=api_key.strip() or None)


def search_transcripts(gene_symbol: str, organism: str = "Homo sapiens", limit: int = 15):
    """Find candidate mRNA RefSeq records for a gene symbol."""
    term = f'({gene_symbol}[Gene]) AND ({organism}[Organism]) AND (biomol_mrna[PROP])'
    with Entrez.esearch(db="nuccore", term=term, retmax=limit, sort="relevance") as handle:
        search_result = Entrez.read(handle)

    ids = search_result.get("IdList", [])
    if not ids:
        return []

    with Entrez.esummary(db="nuccore", id=",".join(ids), retmode="xml") as handle:
        summaries = Entrez.read(handle)

    candidates = []
    for item in summaries:
        accession = str(item.get("AccessionVersion", ""))
        title = str(item.get("Title", ""))
        caption = str(item.get("Caption", ""))
        if accession:
            candidates.append({
                "accession": accession,
                "caption": caption,
                "title": title,
                "label": f"{accession} — {title}"
            })
    return candidates


def fetch_genbank_record(accession: str):
    """Download a GenBank nucleotide record, including its feature annotations."""
    accession = accession.strip()
    if not accession:
        raise ValueError("Enter a RefSeq accession or select a transcript.")
    with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")
    return record


def extract_gene_metadata(record):
    """Collect useful gene and transcript information from GenBank annotations."""
    metadata = {
        "Accession": record.id,
        "Record name": record.name,
        "Description": record.description,
        "Sequence length (nt)": len(record.seq),
        "Organism": record.annotations.get("organism", "Not annotated"),
        "Molecule type": record.annotations.get("molecule_type", "Not annotated"),
        "Topology": record.annotations.get("topology", "Not annotated"),
        "Definition": record.description,
    }

    source = next((f for f in record.features if f.type == "source"), None)
    if source:
        q = source.qualifiers
        metadata["Taxonomy"] = "; ".join(record.annotations.get("taxonomy", []))
        metadata["Source organism"] = ", ".join(q.get("organism", []))

    gene_feature = next((f for f in record.features if f.type == "gene"), None)
    if gene_feature:
        metadata["Gene symbol"] = ", ".join(gene_feature.qualifiers.get("gene", []))

    cds_features = [f for f in record.features if f.type == "CDS"]
    if cds_features:
        cds = cds_features[0]
        metadata["CDS product"] = ", ".join(cds.qualifiers.get("product", [])) or "Not annotated"
        metadata["Protein accession"] = ", ".join(cds.qualifiers.get("protein_id", [])) or "Not annotated"
        metadata["Protein length (aa)"] = len(cds.qualifiers.get("translation", [""])[0]) if cds.qualifiers.get("translation") else "Not provided"
    return metadata


def get_cds_feature(record):
    """Return the first CDS feature and its 1-based inclusive boundaries.

    This MVP supports a contiguous CDS on a transcript record. Complex joined
    CDS locations are rejected rather than silently mapped incorrectly.
    """
    cds_features = [f for f in record.features if f.type == "CDS"]
    if not cds_features:
        raise ValueError(
            "No CDS feature was found in this GenBank record. "
            "Choose a coding RefSeq transcript or provide CDS coordinates manually."
        )

    cds = cds_features[0]
    if getattr(cds.location, "operator", None) == "join":
        raise ValueError(
            "This transcript has a joined/multipart CDS. This MVP does not yet "
            "map multipart CDS coordinates safely."
        )

    start = int(cds.location.start) + 1
    end = int(cds.location.end)
    cds_sequence = str(cds.extract(record.seq)).upper().replace("U", "T")

    return {
        "start": start,
        "end": end,
        "sequence": cds_sequence,
        "feature": cds
    }


def search_variant_database(hgvs: str, max_results: int = 10):
    """Search ClinVar and dbSNP through NCBI Entrez using an HGVS expression.

    A search hit is evidence that a record was found by this query. No hit is
    not proof that a variant is absent from every database.
    """
    results = {}

    for database in ("clinvar", "snp"):
        try:
            with Entrez.esearch(
                db=database,
                term=f'"{hgvs}"',
                retmax=max_results
            ) as handle:
                search_result = Entrez.read(handle)

            ids = search_result.get("IdList", [])
            results[database] = {
                "count": int(search_result.get("Count", 0)),
                "ids": ids,
                "status": "Records found" if ids else "No exact-query hit"
            }
        except Exception as exc:
            results[database] = {
                "count": 0,
                "ids": [],
                "status": f"Search unavailable: {exc}"
            }

    # Retrieve compact ClinVar summaries where possible.
    clinvar_ids = results.get("clinvar", {}).get("ids", [])
    if clinvar_ids:
        try:
            with Entrez.esummary(db="clinvar", id=",".join(clinvar_ids), retmode="xml") as handle:
                summaries = Entrez.read(handle)
            rows = []
            for item in summaries:
                rows.append({
                    "UID": str(item.get("uid", "")),
                    "Title": str(item.get("title", "")),
                    "Accession": str(item.get("accession", "")),
                    "Clinical significance": str(item.get("clinical_significance", {}).get("description", "Not provided"))
                    if isinstance(item.get("clinical_significance", {}), dict)
                    else str(item.get("clinical_significance", "Not provided"))
                })
            results["clinvar"]["summaries"] = rows
        except Exception as exc:
            results["clinvar"]["summary_error"] = str(exc)

    return results
