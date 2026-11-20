import re
import requests
import streamlit as st
import pandas as pd

# ============================================================
# CONFIGURATION
# ============================================================

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
PUBMED_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MYVARIANT_API = "https://myvariant.info/v1"
IPTMNET_API = "https://research.bioinformatics.udel.edu/iptmnet/api/v1"
GTEX_API = "https://gtexportal.org/api/v2"
RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"

REQUEST_TIMEOUT = 25
HEADERS = {
    "User-Agent": "HelixoraBioinformaticsStudio/1.0",
    "Accept": "application/json"
}


# ============================================================
# GENERAL HELPERS
# ============================================================

def api_request(url, params=None, headers=None):
    """GET JSON safely; return a consistent error dictionary on failure."""
    try:
        response = requests.get(
            url,
            params=params,
            headers=headers or HEADERS,
            timeout=REQUEST_TIMEOUT
        )
        response.raise_for_status()
        return response.json()
    except requests.exceptions.HTTPError as error:
        status = error.response.status_code if error.response is not None else "unknown"
        return {"error": f"HTTP {status}: {error}"}
    except requests.exceptions.RequestException as error:
        return {"error": str(error)}
    except ValueError:
        return {"error": "The server returned invalid JSON."}


def safe_text(value, default="Not available"):
    if value is None or value == "":
        return default
    if isinstance(value, (dict, list)):
        return str(value)
    return str(value)


def nested_get(data, *keys, default=None):
    current = data
    for key in keys:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
    return default if current is None else current


def first_value(value):
    if isinstance(value, list):
        return value[0] if value else None
    return value


def get_position(location):
    if not isinstance(location, dict):
        return None
    start = nested_get(location, "start", "value")
    end = nested_get(location, "end", "value")
    if start is None:
        return None
    if end is None or start == end:
        return start
    return f"{start}-{end}"


def to_list(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


# ============================================================
# 1. SNP IDENTIFICATION - MyVariant.info ONLY
# ============================================================

def identify_snp(rsid):
    rsid = rsid.strip().lower()

    if not re.fullmatch(r"rs\d+", rsid):
        return {"error": "Invalid SNP ID. Example: rs1260326"}

    data = api_request(
        f"{MYVARIANT_API}/query",
        params={
            "q": f"dbsnp.rsid:{rsid}",
            "size": 20
        }
    )

    if "error" in data:
        return data

    hits = data.get("hits", [])
    exact_hits = []

    for hit in hits:
        returned_rsid = nested_get(hit, "dbsnp", "rsid", default="")
        if str(returned_rsid).lower() == rsid:
            exact_hits.append(hit)

    if not exact_hits:
        return {
            "error": (
                f"{rsid} was not found as an exact dbSNP match in MyVariant.info. "
                "Try again later or verify the rsID."
            )
        }

    # MyVariant can return multiple genomic records for one rsID.
    # Prefer a record containing GCKR-related annotations when available.
    hit = next(
        (
            item for item in exact_hits
            if "GCKR" in str(nested_get(item, "dbnsfp", "genename", default="")).upper()
            or "GCKR" in str(nested_get(item, "cadd", "gene", "genename", default="")).upper()
        ),
        exact_hits[0]
    )

    dbsnp = hit.get("dbsnp", {}) or {}
    dbnsfp = hit.get("dbnsfp", {}) or {}
    cadd = hit.get("cadd", {}) or {}

    ref_allele = safe_text(dbsnp.get("ref"))
    alt_value = dbsnp.get("alt")
    if isinstance(alt_value, list):
        alt_allele = ", ".join(str(x) for x in alt_value)
    else:
        alt_allele = safe_text(alt_value)

    gene_name = dbnsfp.get("genename", "Not available")
    if isinstance(gene_name, list):
        gene_name = ", ".join(map(str, gene_name))

    ensembl_annotation = first_value(dbnsfp.get("ensembl", {})) or {}
    aa_annotation = first_value(dbnsfp.get("aa", {})) or {}

    protein_change = "Not available"
    if isinstance(aa_annotation, dict):
        aa_ref = aa_annotation.get("ref")
        aa_alt = aa_annotation.get("alt")
        aa_pos = aa_annotation.get("pos")
        if aa_ref and aa_alt and aa_pos:
            protein_change = f"{aa_ref}{aa_pos}{aa_alt}"

    # Genomic mappings may be represented differently between MyVariant records.
    genomic_results = []
    for assembly in ("hg19", "hg38"):
        loc = dbsnp.get(assembly) or hit.get(assembly)
        if isinstance(loc, dict):
            genomic_results.append({
                "Assembly": assembly.upper(),
                "Chromosome": loc.get("chrom", dbsnp.get("chrom", hit.get("chrom", "Unknown"))),
                "Start": loc.get("start", "Unknown"),
                "End": loc.get("end", "Unknown"),
                "Reference": ref_allele,
                "Alternate": alt_allele
            })

    # Extract consequence hints from MyVariant annotations.
    # These are not equivalent to a full Ensembl VEP transcript consequence report.
    consequence_rows = []
    for source_name, source in (("dbNSFP", dbnsfp), ("CADD", cadd)):
        if not isinstance(source, dict):
            continue

        gene = source.get("genename")
        if isinstance(gene, list):
            gene = first_value(gene)

        gene_id = nested_get(source, "ensembl", "geneid")
        transcript_id = nested_get(source, "ensembl", "transcriptid")

        if isinstance(source.get("ensembl"), list):
            ens = first_value(source.get("ensembl")) or {}
            gene_id = gene_id or ens.get("geneid")
            transcript_id = transcript_id or ens.get("transcriptid")

        consequence = (
            source.get("consequence")
            or source.get("annotype")
            or source.get("effect")
            or "Not available"
        )

        if gene or gene_id or transcript_id or consequence != "Not available":
            consequence_rows.append({
                "Annotation Source": source_name,
                "Gene": safe_text(gene),
                "Ensembl Gene ID": safe_text(gene_id),
                "Transcript ID": safe_text(transcript_id),
                "Consequence": safe_text(consequence),
                "Amino Acid Change": protein_change
                if source_name == "dbNSFP" else "Not available"
            })

    gckr_consequences = [
        row for row in consequence_rows
        if row["Gene"].upper() == "GCKR"
    ]
    other_consequences = [
        row for row in consequence_rows
        if row["Gene"].upper() not in ("GCKR", "NOT AVAILABLE")
    ]

    clinvar = hit.get("clinvar", {}) or {}
    clinvar_records = []
    if isinstance(clinvar, dict):
        clinvar_records = to_list(clinvar.get("rcv", []))

    # Extract allele frequency annotations when available.
    frequency_records = []
    for key in ("gnomad_exome", "gnomad_genome", "exac", " allele"):
        value = dbsnp.get(key) or dbnsfp.get(key)
        if value:
            frequency_records.append({
                "Database": key.strip(),
                "Data": safe_text(value)
            })

    return {
        "SNP ID": rsid,
        "Gene": safe_text(gene_name),
        "Ensembl Gene ID": safe_text(ensembl_annotation.get("geneid")),
        "Reference Allele": ref_allele,
        "Alternate Allele": alt_allele,
        "Protein Change": protein_change,
        "Most Severe Consequence": safe_text(
            dbnsfp.get("consequence") or cadd.get("annotype")
        ),
        "Genomic Mapping": genomic_results,
        "GCKR Transcript Consequences": gckr_consequences,
        "Other Transcript Consequences": other_consequences,
        "Regulatory Consequences": [],
        "Allele Frequency Annotations": frequency_records,
        "ClinVar Information": clinvar_records,
        "dbSNP": f"https://www.ncbi.nlm.nih.gov/snp/{rsid}",
        "ClinVar": f"https://www.ncbi.nlm.nih.gov/clinvar/?term={rsid}",
        "Ensembl": (
            "https://www.ensembl.org/Homo_sapiens/Variation/Explore"
            f"?v={rsid}"
        ),
        "Annotation Source": "MyVariant.info",
        "VEP Status": (
            "Not queried. Ensembl REST is intentionally disabled because "
            "the user's network receives HTTP 403."
        )
    }


# ============================================================
# 2. GENE AND PROTEIN INFORMATION - UniProt
# ============================================================

def get_gene_information(gene_symbol):
    data = api_request(
        f"{UNIPROT_API}/search",
        params={
            "query": f"gene_exact:{gene_symbol} AND organism_id:9606",
            "format": "json",
            "size": 20
        }
    )

    if "error" in data:
        return data

    results = []
    for entry in data.get("results", []):
        genes = entry.get("genes", [])
        names = [
            nested_get(item, "geneName", "value", default="")
            for item in genes
        ]

        if gene_symbol.upper() not in [name.upper() for name in names if name]:
            continue

        protein_name = nested_get(
            entry, "proteinDescription", "recommendedName",
            "fullName", "value", default="Unknown"
        )

        results.append({
            "Protein name": protein_name,
            "Gene": ", ".join(name for name in names if name),
            "UniProt ID": entry.get("primaryAccession"),
            "Length": nested_get(entry, "sequence", "length"),
            "Organism": nested_get(
                entry, "organism", "scientificName", default="Unknown"
            )
        })

    return results


def get_protein_information(uniprot_id):
    if not re.fullmatch(r"[A-Z0-9]{6,10}", uniprot_id.upper()):
        return {"error": "Invalid UniProt accession."}
    return api_request(f"{UNIPROT_API}/{uniprot_id.upper()}.json")


# ============================================================
# 3. PTMs: UniProt + iPTMnet
# ============================================================

PTM_FEATURE_TYPES = {
    "MODIFIED RESIDUE", "LIPIDATION", "GLYCOSYLATION",
    "CROSSLINK", "DISULFID"
}


def _normalise_evidence(evidences):
    """Return a readable evidence count without assuming one JSON shape."""
    if evidences is None:
        return 0
    if isinstance(evidences, list):
        return len(evidences)
    if isinstance(evidences, dict):
        return len(evidences)
    return 1


def _sequence_residue(protein_data, position):
    """Get the amino acid at a 1-based position from the UniProt sequence."""
    sequence = nested_get(protein_data, "sequence", "value", default="")
    try:
        pos = int(position)
        if sequence and 1 <= pos <= len(sequence):
            return sequence[pos - 1]
    except (TypeError, ValueError):
        pass
    return "Not available"


def get_uniprot_ptms(protein_data):
    """Extract experimentally/curator annotated PTM-like UniProt features."""
    if not isinstance(protein_data, dict) or "error" in protein_data:
        return []

    results = []
    for feature in protein_data.get("features", []):
        feature_type = str(feature.get("type", "")).upper()
        if feature_type not in PTM_FEATURE_TYPES:
            continue

        start = nested_get(feature, "location", "start", "value")
        end = nested_get(feature, "location", "end", "value")
        position = get_position(feature.get("location", {}))
        residue = _sequence_residue(protein_data, start) if start is not None else "Not available"
        description = feature.get("description") or "Not specified"
        results.append({
            "Source": "UniProt",
            "Modification": feature_type.title(),
            "Position": position if position is not None else "Not available",
            "Residue": residue,
            "Description": description,
            "Enzyme": "Not specified",
            "Evidence count": _normalise_evidence(feature.get("evidences")),
            "Evidence / PMID": "See UniProt feature evidence",
            "Evidence class": "UniProt annotation; inspect evidence codes",
            "Score": "Not provided"
        })
    return results


def _extract_iptmnet_table(payload):
    """iPTMnet v1 substrate endpoint returns {form, table:[...]}.

    Keep a defensive fallback for alternate wrappers, but parse the documented
    table schema first rather than recursively mistaking metadata for sites.
    """
    if isinstance(payload, dict):
        table = payload.get("table")
        if isinstance(table, list):
            return table
        # Some deployments wrap the documented response.
        for key in ("data", "results", "result"):
            value = payload.get(key)
            if isinstance(value, dict) and isinstance(value.get("table"), list):
                return value["table"]
            if isinstance(value, list):
                return value
    if isinstance(payload, list):
        return payload
    return []


def _join_values(value):
    if value is None:
        return "Not available"
    if isinstance(value, list):
        return "; ".join(_join_values(item) for item in value) if value else "Not available"
    if isinstance(value, dict):
        # Evidence/source objects often contain a name and/or URL.
        name = value.get("name") or value.get("id") or value.get("label")
        if name:
            return str(name)
        return "; ".join(f"{k}: {v}" for k, v in value.items())
    return str(value)


def get_iptmnet_ptms(uniprot_id):
    """Retrieve documented PTM substrate sites from the iPTMnet v1 API."""
    if not re.fullmatch(r"[A-Z0-9]{6,10}", str(uniprot_id).upper()):
        return [], "Invalid UniProt accession supplied to iPTMnet."

    # The documented API route is /v1/{id}/substrate/.
    payload = api_request(
        f"{IPTMNET_API}/{uniprot_id.upper()}/substrate/"
    )
    if isinstance(payload, dict) and "error" in payload:
        return [], payload["error"]

    raw_records = _extract_iptmnet_table(payload)
    results = []
    seen = set()

    for item in raw_records:
        if not isinstance(item, dict):
            continue
        site = item.get("site") or item.get("site_position") or "Not available"
        ptm_type = item.get("PTM type") or item.get("ptm_type") or item.get("ptm") or item.get("modification") or "Not available"
        enzyme = item.get("enzyme") or item.get("ptm_enzyme") or {}
        if isinstance(enzyme, dict):
            enzyme_text = enzyme.get("name") or enzyme.get("id") or enzyme.get("type") or "Not specified"
        else:
            enzyme_text = _join_values(enzyme)

        sources = item.get("sources") or item.get("source") or "Not available"
        pmids = item.get("PMIDs") or item.get("pmids") or item.get("pmid") or []
        score = item.get("score", "Not provided")
        position = item.get("site_position")
        if position is None:
            match = re.search(r"\d+", str(site))
            position = match.group(0) if match else "Not available"
        residue_match = re.match(r"([A-Z])", str(site))
        residue = residue_match.group(1) if residue_match else "Not available"

        row = {
            "Source": "iPTMnet",
            "Modification": str(ptm_type),
            "Position": str(position),
            "Residue": residue,
            "Description": str(site),
            "Enzyme": str(enzyme_text),
            "Evidence count": len(to_list(sources)) if sources not in (None, "Not available") else 0,
            "Evidence / PMID": f"Sources: {_join_values(sources)} | PMID: {_join_values(pmids)}",
            "Evidence class": "Integrated iPTMnet record; may combine curated and text-mined evidence",
            "Score": score
        }
        signature = tuple(str(row[k]) for k in ("Modification", "Position", "Residue", "Enzyme", "Evidence / PMID"))
        if signature not in seen:
            seen.add(signature)
            results.append(row)

    return results, None


def get_ptms(protein_data, uniprot_id="Q14397"):
    """Merge UniProt and iPTMnet records while retaining source provenance."""
    uniprot_records = get_uniprot_ptms(protein_data)
    iptm_records, iptm_error = get_iptmnet_ptms(uniprot_id)
    combined = uniprot_records + iptm_records

    # Avoid exact duplicate records from the same source.
    unique = []
    seen = set()
    for row in combined:
        signature = (row.get("Source"), row.get("Modification"), row.get("Position"), row.get("Residue"), row.get("Enzyme"))
        if signature not in seen:
            seen.add(signature)
            unique.append(row)
    return unique, iptm_error


# ============================================================
# 4. PROTEIN DOMAINS: UniProt + InterPro
# ============================================================

def get_uniprot_domains(protein_data):
    """Extract UniProt DOMAIN/REGION/REPEAT features with sequence coordinates."""
    if not isinstance(protein_data, dict) or "error" in protein_data:
        return []

    results = []
    for feature in protein_data.get("features", []):
        kind = str(feature.get("type", "")).upper()
        if kind not in ("DOMAIN", "REGION", "REPEAT"):
            continue
        start = nested_get(feature, "location", "start", "value")
        end = nested_get(feature, "location", "end", "value")
        try:
            length = int(end) - int(start) + 1
        except (TypeError, ValueError):
            length = None
        results.append({
            "Name": feature.get("description") or kind.title(),
            "Type": kind.title(),
            "Accession": "UniProt feature",
            "Start": start,
            "End": end,
            "Length (aa)": length,
            "Database": "UniProt",
            "Integrated entry": "Not applicable",
            "Evidence": "UniProt feature annotation"
        })
    return results


def get_interpro_domains(uniprot_id, max_pages=5):
    """Fetch InterPro matches for a UniProt protein, following API pagination."""
    if not re.fullmatch(r"[A-Z0-9]{6,10}", str(uniprot_id).upper()):
        return [], "Invalid UniProt accession supplied to InterPro."

    url = f"{INTERPRO_API}/entry/interpro/protein/uniprot/{uniprot_id.upper()}/"
    params = {"page_size": 100}
    records = []
    next_url = url
    first_error = None

    for _ in range(max_pages):
        payload = api_request(next_url, params=params if next_url == url else None,
                              headers={**HEADERS, "Accept": "application/json"})
        params = None
        if isinstance(payload, dict) and "error" in payload:
            first_error = payload["error"]
            break
        page_records = payload.get("results", []) if isinstance(payload, dict) else []
        records.extend(page_records)
        next_url = payload.get("next") if isinstance(payload, dict) else None
        if not next_url:
            break

    parsed = []
    for item in records:
        metadata = item.get("metadata", {}) if isinstance(item, dict) else {}
        accession = metadata.get("accession", "Not available")
        name = metadata.get("name", "Not available")
        entry_type = metadata.get("type", "Not specified")
        source_db = metadata.get("source_database", "InterPro")
        integrated = metadata.get("integrated", "")

        # API records commonly carry protein locations under proteins[].
        locations = []
        for protein in item.get("proteins", []) if isinstance(item, dict) else []:
            for entry_location in protein.get("entry_protein_locations", []) or []:
                for fragment in entry_location.get("fragments", []) or []:
                    locations.append((fragment.get("start"), fragment.get("end")))
        # Some API payload versions place fragments directly on the result.
        if not locations:
            for entry_location in item.get("entry_protein_locations", []) if isinstance(item, dict) else []:
                for fragment in entry_location.get("fragments", []) or []:
                    locations.append((fragment.get("start"), fragment.get("end")))

        if not locations:
            locations = [(None, None)]
        for start, end in locations:
            try:
                length = int(end) - int(start) + 1
            except (TypeError, ValueError):
                length = None
            parsed.append({
                "Name": name,
                "Type": str(entry_type).title(),
                "Accession": accession,
                "Start": start,
                "End": end,
                "Length (aa)": length,
                "Database": source_db,
                "Integrated entry": integrated or "Not integrated",
                "Evidence": "InterPro protein match"
            })

    # Deduplicate identical entry/location rows.
    unique = []
    seen = set()
    for row in parsed:
        sig = (row["Accession"], row["Start"], row["End"], row["Database"])
        if sig not in seen:
            seen.add(sig)
            unique.append(row)
    return unique, first_error


def get_protein_domains(protein_data, uniprot_id="Q14397"):
    """Combine UniProt features with InterPro matches; return provenance and errors."""
    uniprot_rows = get_uniprot_domains(protein_data)
    interpro_rows, interpro_error = get_interpro_domains(uniprot_id)
    return uniprot_rows + interpro_rows, interpro_error


def get_protein_family(protein_data):
    if not isinstance(protein_data, dict) or "error" in protein_data:
        return []

    results = []
    for ref in protein_data.get("uniProtKBCrossReferences", []):
        database = ref.get("database")
        if database in ("InterPro", "Pfam", "PANTHER", "SUPFAM"):
            props = ref.get("properties", [])
            description = next(
                (p.get("value") for p in props if p.get("key") in ("EntryName", "Description")),
                "Not specified"
            )
            results.append({"Database": database, "ID": ref.get("id"), "Description": description})
    return results

# ============================================================
# 5. TISSUE EXPRESSION - GTEx API V2, NO ENSEMBL REST
# ============================================================

def get_gtex_gene_record(gene_symbol):
    """
    Resolve a gene symbol directly through GTEx reference/geneSearch.
    No Ensembl REST request is made.
    """
    data = api_request(
        f"{GTEX_API}/reference/geneSearch",
        params={
            "geneId": gene_symbol,
            "gencodeVersion": "v39",
            "genomeBuild": "GRCh38/hg38",
            "page": 0,
            "itemsPerPage": 50
        }
    )

    if "error" in data:
        return None, data["error"]

    records = data.get("data", []) if isinstance(data, dict) else []
    if not records:
        return None, f"GTEx returned no gene reference record for {gene_symbol}."

    exact = next(
        (
            item for item in records
            if str(item.get("geneSymbol", "")).upper() == gene_symbol.upper()
        ),
        None
    )

    if exact is None:
        return None, (
            f"GTEx did not return an exact gene-symbol match for {gene_symbol}. "
            "Returned records were not used to avoid selecting the wrong gene."
        )

    gencode_id = (
        exact.get("gencodeId")
        or exact.get("gencodeIdVersion")
        or exact.get("geneId")
    )

    if not gencode_id:
        return None, "GTEx gene record did not contain a GENCODE identifier."

    # GTEx expression endpoints expect versioned GENCODE IDs.
    if re.fullmatch(r"ENSG\d+", str(gencode_id)):
        return None, (
            f"GTEx returned an unversioned gene ID ({gencode_id}); "
            "a versioned GENCODE ID is required for this expression query."
        )

    return {
        "gene_symbol": gene_symbol,
        "gencode_id": gencode_id,
        "gene_record": exact
    }, None


def get_gtex_expression(gene_symbol):
    """
    Retrieve GTEx median gene expression across tissues.
    Uses GTEx API v2 reference/geneSearch and expression/medianGeneExpression.
    """
    gene_record, error = get_gtex_gene_record(gene_symbol)
    if error:
        return [], error

    gencode_id = gene_record["gencode_id"]

    data = api_request(
        f"{GTEX_API}/expression/medianGeneExpression",
        params={
            "gencodeId": gencode_id,
            "datasetId": "gtex_v10",
            "page": 0,
            "itemsPerPage": 250
        }
    )

    if "error" in data:
        return [], data["error"]

    records = data.get("data", []) if isinstance(data, dict) else []
    results = []

    for item in records:
        tissue = (
            item.get("tissueSiteDetail")
            or item.get("tissueSiteDetailId")
            or item.get("tissue")
            or "Unknown"
        )

        median = (
            item.get("median")
            if item.get("median") is not None
            else item.get("medianExpression")
        )

        results.append({
            "Tissue": tissue,
            "Median TPM": median if median is not None else "Not available",
            "Gene": gene_symbol,
            "GENCODE ID": gencode_id,
            "Dataset": "GTEx v10"
        })

    results.sort(
        key=lambda row: (
            row["Median TPM"]
            if isinstance(row["Median TPM"], (int, float))
            else -1
        ),
        reverse=True
    )

    if not results:
        return [], (
            "GTEx resolved the gene, but no median-expression records were "
            "returned for this GENCODE ID and dataset."
        )

    return results, None


def get_tissue_specificity(protein_data):
    results = []
    if not isinstance(protein_data, dict) or "error" in protein_data:
        return results

    for comment in protein_data.get("comments", []):
        if comment.get("commentType") == "TISSUE SPECIFICITY":
            for item in comment.get("texts", []):
                results.append(item.get("value", ""))
    return results


# ============================================================
# 6. PUBMED LITERATURE SEARCH
# ============================================================

def search_pubmed(gene_symbol, limit=10):
    term = (
        f'("{gene_symbol}"[Title/Abstract] OR '
        f'"glucokinase regulatory protein"[Title/Abstract])'
    )

    search = api_request(
        f"{PUBMED_API}/esearch.fcgi",
        params={
            "db": "pubmed",
            "term": term,
            "retmode": "json",
            "retmax": limit,
            "sort": "relevance"
        }
    )

    if "error" in search:
        return []

    ids = nested_get(search, "esearchresult", "idlist", default=[])
    if not ids:
        return []

    summary = api_request(
        f"{PUBMED_API}/esummary.fcgi",
        params={
            "db": "pubmed",
            "id": ",".join(ids),
            "retmode": "json"
        }
    )

    if "error" in summary:
        return []

    results = []
    for pmid in ids:
        article = nested_get(summary, "result", pmid, default={})
        results.append({
            "Title": article.get("title", "Unknown"),
            "Authors": ", ".join(
                author.get("name", "")
                for author in article.get("authors", [])[:5]
            ),
            "Journal": article.get("fulljournalname", "Unknown"),
            "Year": article.get("pubdate", "Unknown"),
            "PMID": pmid,
            "PubMed": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        })

    return results


# ============================================================
# 7. PDB STRUCTURES
# ============================================================

def search_pdb(uniprot_id):
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_accession"
                ),
                "operator": "exact_match",
                "value": uniprot_id
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 30}
        }
    }

    try:
        response = requests.post(
            RCSB_SEARCH_API,
            json=query,
            headers=HEADERS,
            timeout=REQUEST_TIMEOUT
        )
        response.raise_for_status()
        return [
            item["identifier"]
            for item in response.json().get("result_set", [])
            if "identifier" in item
        ]
    except (requests.exceptions.RequestException, ValueError, KeyError):
        return []


def download_pdb(pdb_id):
    if not re.fullmatch(r"[A-Za-z0-9]{4}", pdb_id):
        return None

    try:
        response = requests.get(
            f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb",
            timeout=REQUEST_TIMEOUT
        )
        response.raise_for_status()
        return response.content
    except requests.exceptions.RequestException:
        return None


# ============================================================
# MAIN STREAMLIT PAGE
# ============================================================

def render_gckr_page():
    st.title("Gene & Protein Characterization")
    st.write(
        "Explore SNPs, gene and protein annotations, PTMs, domains, "
        "protein family, tissue expression, literature and 3D structures."
    )
    st.caption(
        "Sources: MyVariant.info, Ensembl Website, UniProt, "
        "iPTMnet, GTEx, PubMed and RCSB PDB."
    )
    st.info(
        "Ensembl REST API is not used in this version because requests from "
        "the current network returned HTTP 403. Ensembl is provided as a "
        "manual external record link."
    )
    st.divider()

    # 1. SNP Identification
    st.header("1. SNP Identification")
    rsid = st.text_input(
        "Enter SNP ID",
        value="rs1260326",
        key="gckr_rsid"
    )

    if st.button("Identify SNP", key="identify_snp"):
        with st.spinner("Retrieving SNP annotations from MyVariant.info..."):
            st.session_state["snp_data"] = identify_snp(rsid)

    if "snp_data" in st.session_state:
        snp = st.session_state["snp_data"]

        if "error" in snp:
            st.error(snp["error"])
        else:
            st.success(f"Retrieved {snp.get('SNP ID', rsid)}")

            c1, c2, c3 = st.columns(3)
            c1.metric("SNP ID", snp.get("SNP ID", "Unknown"))
            c2.metric(
                "Reference Allele",
                snp.get("Reference Allele", "Not available")
            )
            c3.metric(
                "Alternate Allele",
                snp.get("Alternate Allele", "Not available")
            )

            st.write("**Gene annotation:**", snp.get("Gene", "Not available"))
            st.write(
                "**Ensembl Gene ID (if supplied by MyVariant):**",
                snp.get("Ensembl Gene ID", "Not available")
            )
            st.write(
                "**Most severe consequence annotation:**",
                snp.get("Most Severe Consequence", "Not available")
            )
            st.write(
                "**Protein change annotation:**",
                snp.get("Protein Change", "Not available")
            )
            st.caption(
                "These consequence fields are annotations returned by "
                "MyVariant.info. They are not a full Ensembl VEP report."
            )
            st.caption(
                "Reference and alternate alleles are not the same as "
                "ancestral and minor alleles. Minor allele frequency depends "
                "on population."
            )

            st.subheader("Genomic Coordinates")
            if snp.get("Genomic Mapping"):
                st.dataframe(
                    pd.DataFrame(snp["Genomic Mapping"]),
                    use_container_width=True,
                    hide_index=True
                )
            else:
                st.info(
                    "No hg19/hg38 genomic mapping was returned in the "
                    "selected MyVariant record."
                )

            st.subheader("GCKR Consequence Annotations")
            if snp.get("GCKR Transcript Consequences"):
                st.dataframe(
                    pd.DataFrame(snp["GCKR Transcript Consequences"]),
                    use_container_width=True,
                    hide_index=True
                )
            else:
                st.info(
                    "MyVariant.info did not return a GCKR-specific transcript "
                    "consequence row for this record. This is not evidence "
                    "that the variant has no transcript consequence."
                )

            with st.expander("Other consequence annotations"):
                other = snp.get("Other Transcript Consequences", [])
                if other:
                    st.dataframe(
                        pd.DataFrame(other),
                        use_container_width=True,
                        hide_index=True
                    )
                else:
                    st.write("No additional consequence annotations returned.")

            with st.expander("ClinVar records returned by MyVariant.info"):
                records = snp.get("ClinVar Information", [])
                if records:
                    st.json(records)
                else:
                    st.write(
                        "No ClinVar record was included in this MyVariant "
                        "response. Use the external ClinVar search."
                    )

            with st.expander("Allele frequency annotations"):
                frequency = snp.get("Allele Frequency Annotations", [])
                if frequency:
                    st.dataframe(
                        pd.DataFrame(frequency),
                        use_container_width=True,
                        hide_index=True
                    )
                else:
                    st.write(
                        "No allele-frequency field was found in the selected "
                        "MyVariant record. Do not infer minor allele status."
                    )

            st.write("**Annotation source:**", snp.get("Annotation Source"))
            st.caption(snp.get("VEP Status", ""))

            st.subheader("External Database Records")
            st.markdown(f"[Open dbSNP record]({snp['dbSNP']})")
            st.markdown(f"[Search ClinVar]({snp['ClinVar']})")
            st.markdown(f"[Open Ensembl record]({snp['Ensembl']})")

    st.divider()

    # 2. Gene and protein
    st.header("2. Gene & Protein Information")
    gene_symbol = st.text_input(
        "Gene symbol",
        value="GCKR",
        key="gckr_gene"
    )
    uniprot_id = st.text_input(
        "UniProt accession",
        value="Q14397",
        key="gckr_uniprot"
    )

    if st.button(
        "Retrieve Gene and Protein Information",
        key="retrieve_gene_protein"
    ):
        with st.spinner("Retrieving UniProt information..."):
            st.session_state["gene_data"] = get_gene_information(gene_symbol)
            st.session_state["protein_data"] = get_protein_information(uniprot_id)

    protein_data = st.session_state.get("protein_data")

    if "gene_data" in st.session_state:
        gene_data = st.session_state["gene_data"]
        st.subheader("Gene-associated protein records")

        if isinstance(gene_data, dict) and "error" in gene_data:
            st.error(gene_data["error"])
        elif gene_data:
            st.dataframe(
                pd.DataFrame(gene_data),
                use_container_width=True,
                hide_index=True
            )
        else:
            st.warning("No matching human UniProt protein records found.")

    if protein_data is not None:
        if "error" in protein_data:
            st.error(protein_data["error"])
        else:
            protein_name = nested_get(
                protein_data, "proteinDescription", "recommendedName",
                "fullName", "value", default="Unknown"
            )
            length = nested_get(
                protein_data, "sequence", "length", default="Unknown"
            )

            c1, c2 = st.columns(2)
            c1.metric("Protein Length (aa)", length)
            c2.metric(
                "UniProt Accession",
                protein_data.get("primaryAccession", "Unknown")
            )

            st.write("**Protein name:**", protein_name)
            st.write(
                "**Organism:**",
                nested_get(
                    protein_data, "organism", "scientificName",
                    default="Unknown"
                )
            )
            st.markdown(
                f"[Open UniProt entry]"
                f"(https://www.uniprot.org/uniprotkb/{uniprot_id})"
            )

    st.divider()

    # 3. PTM
    st.header("3. Post-translational Modifications (PTMs)")
    st.write(
        "PTM sites are collected independently from UniProt and iPTMnet. "
        "The source of each annotation is retained so evidence is not conflated."
    )
    if isinstance(protein_data, dict) and "error" not in protein_data:
        if st.button("Retrieve and analyze PTMs", key="retrieve_ptms"):
            with st.spinner("Querying UniProt and iPTMnet..."):
                ptm_rows, ptm_error = get_ptms(protein_data, uniprot_id)
                st.session_state["ptm_data"] = ptm_rows
                st.session_state["ptm_error"] = ptm_error

        if "ptm_data" in st.session_state:
            ptms = st.session_state["ptm_data"]
            ptm_error = st.session_state.get("ptm_error")
            if ptm_error:
                st.warning(f"iPTMnet could not be queried: {ptm_error}")
            if ptms:
                ptm_df = pd.DataFrame(ptms)
                c1, c2, c3 = st.columns(3)
                c1.metric("Total annotated records", len(ptm_df))
                c2.metric("UniProt records", int((ptm_df["Source"] == "UniProt").sum()))
                c3.metric("iPTMnet records", int((ptm_df["Source"] == "iPTMnet").sum()))
                st.subheader("PTM annotation table")
                st.dataframe(ptm_df, use_container_width=True, hide_index=True)
                st.download_button(
                    "Download PTM annotations (CSV)",
                    data=ptm_df.to_csv(index=False).encode("utf-8-sig"),
                    file_name=f"{uniprot_id}_PTM_annotations.csv",
                    mime="text/csv",
                    key="download_ptm_csv"
                )
                st.subheader("PTM distribution by modification type")
                counts = ptm_df.groupby(["Modification", "Source"]).size().reset_index(name="Records")
                st.bar_chart(counts, x="Modification", y="Records", color="Source")
                numeric_positions = pd.to_numeric(ptm_df["Position"], errors="coerce")
                plot_df = ptm_df.loc[numeric_positions.notna(), ["Modification", "Source"]].copy()
                plot_df["Position"] = numeric_positions[numeric_positions.notna()].astype(int).values
                if not plot_df.empty:
                    st.subheader("PTM site positions")
                    st.scatter_chart(plot_df, x="Position", y="Modification", color="Source")
            else:
                st.info("No PTM sites were returned. This does not prove that GKRP has no PTMs.")
            st.caption(
                "UniProt and iPTMnet annotations are database-level evidence, "
                "not proof that each modification occurs in every tissue or condition. "
                "iPTMnet may integrate curated and text-mined sources; inspect its source/PMID fields."
            )
            st.markdown("[iPTMnet API documentation](https://research.bioinformatics.udel.edu/iptmnet/api/doc/)")
    else:
        st.info("Retrieve valid protein information above before querying PTMs.")

    st.divider()

    # 4. Domains
    st.header("4. Protein Domains")
    st.write(
        "Domain annotations combine UniProt sequence features with InterPro "
        "protein matches. The source database and sequence coordinates are retained."
    )
    if isinstance(protein_data, dict) and "error" not in protein_data:
        if st.button("Retrieve and map protein domains", key="retrieve_domains"):
            with st.spinner("Retrieving UniProt and InterPro domain annotations..."):
                domain_rows, domain_error = get_protein_domains(protein_data, uniprot_id)
                st.session_state["domain_data"] = domain_rows
                st.session_state["domain_error"] = domain_error

        if "domain_data" in st.session_state:
            domains = st.session_state["domain_data"]
            domain_error = st.session_state.get("domain_error")
            if domain_error:
                st.warning(f"InterPro could not be queried: {domain_error}")
            if domains:
                domain_df = pd.DataFrame(domains)
                c1, c2 = st.columns(2)
                c1.metric("Total domain/family annotations", len(domain_df))
                c2.metric("Records with coordinates", int(domain_df["Start"].notna().sum()))
                st.subheader("Domain annotation table")
                st.dataframe(domain_df, use_container_width=True, hide_index=True)
                st.download_button(
                    "Download domain annotations (CSV)",
                    data=domain_df.to_csv(index=False).encode("utf-8-sig"),
                    file_name=f"{uniprot_id}_protein_domains.csv",
                    mime="text/csv",
                    key="download_domain_csv"
                )
                # A Vega-Lite interval map is native to Streamlit and avoids an extra plotting dependency.
                map_df = domain_df.copy()
                map_df["Start_num"] = pd.to_numeric(map_df["Start"], errors="coerce")
                map_df["End_num"] = pd.to_numeric(map_df["End"], errors="coerce")
                map_df = map_df.dropna(subset=["Start_num", "End_num"])
                if not map_df.empty:
                    map_df["Feature"] = map_df["Name"].astype(str) + " | " + map_df["Accession"].astype(str)
                    st.subheader("Protein domain map")
                    st.vega_lite_chart(
                        map_df,
                        {
                            "mark": {"type": "bar", "height": 18},
                            "encoding": {
                                "y": {"field": "Feature", "type": "nominal", "sort": "-x", "title": "Domain / feature"},
                                "x": {"field": "Start_num", "type": "quantitative", "title": "Amino-acid position", "scale": {"zero": False}},
                                "x2": {"field": "End_num"},
                                "color": {"field": "Database", "type": "nominal", "title": "Source"},
                                "tooltip": [
                                    {"field": "Name", "type": "nominal"},
                                    {"field": "Accession", "type": "nominal"},
                                    {"field": "Start_num", "type": "quantitative", "title": "Start"},
                                    {"field": "End_num", "type": "quantitative", "title": "End"},
                                    {"field": "Length (aa)", "type": "quantitative"},
                                    {"field": "Database", "type": "nominal"}
                                ]
                            }
                        },
                        use_container_width=True
                    )
                else:
                    st.info("Domain records were found, but no usable start/end coordinates were supplied for the map.")
            else:
                st.info("No domain or family matches were returned by UniProt/InterPro.")
            st.caption(
                "InterPro matches can include domains, families, repeats and other signature types; "
                "check the Type and Database columns before interpreting a match as a structural domain."
            )
            st.markdown("[InterPro API documentation](https://interpro-documentation.readthedocs.io/en/latest/api.html)")

    st.divider()

    # 5. Family
    st.header("5. Protein Family & Cross-references")
    if isinstance(protein_data, dict) and "error" not in protein_data:
        families = get_protein_family(protein_data)
        if families:
            st.dataframe(
                pd.DataFrame(families),
                use_container_width=True,
                hide_index=True
            )
        else:
            st.info(
                "No InterPro/Pfam/PANTHER/SUPFAM cross-references were found "
                "in this UniProt record."
            )

    st.divider()

    # 6. Tissue expression
    st.header("6. Tissue Expression")
    st.write(
        "GTEx provides quantitative gene-expression measurements across "
        "human tissues. This version resolves the gene directly through "
        "GTEx API v2 and does not call Ensembl REST."
    )

    if st.button(
        "Retrieve GTEx tissue expression",
        key="retrieve_gtex"
    ):
        with st.spinner("Resolving GCKR in GTEx and retrieving expression..."):
            expression, error = get_gtex_expression(gene_symbol)
            st.session_state["gtex_expression"] = expression
            st.session_state["gtex_error"] = error

    if "gtex_expression" in st.session_state:
        expression = st.session_state["gtex_expression"]
        error = st.session_state.get("gtex_error")

        if error:
            st.warning(error)
        elif expression:
            st.dataframe(
                pd.DataFrame(expression),
                use_container_width=True,
                hide_index=True
            )
            st.caption(
                "GTEx v10 median gene-expression values. Tissues are sorted "
                "from highest to lowest when numeric values are returned."
            )
        else:
            st.info("No GTEx expression records were returned.")

    if isinstance(protein_data, dict) and "error" not in protein_data:
        tissue_comments = get_tissue_specificity(protein_data)
        with st.expander("UniProt qualitative tissue-specificity annotation"):
            if tissue_comments:
                for comment in tissue_comments:
                    st.write("•", comment)
            else:
                st.write(
                    "No tissue-specificity comment in this UniProt record."
                )

    st.divider()

    # 7. Literature
    st.header("7. Literature Review")
    if st.button("Search PubMed", key="search_pubmed"):
        with st.spinner("Searching PubMed..."):
            st.session_state["articles"] = search_pubmed(
                gene_symbol,
                limit=15
            )

    if "articles" in st.session_state:
        articles = st.session_state["articles"]
        if articles:
            st.dataframe(
                pd.DataFrame(articles),
                use_container_width=True,
                hide_index=True
            )
        else:
            st.warning("No PubMed articles were returned for this query.")

    st.divider()

    # 8. PDB
    st.header("8. 3D Protein Structure")
    if st.button("Search PDB Structures", key="search_pdb"):
        with st.spinner(
            "Searching RCSB PDB for UniProt-linked structures..."
        ):
            st.session_state["pdb_ids"] = search_pdb(uniprot_id)

    if "pdb_ids" in st.session_state:
        pdb_ids = st.session_state["pdb_ids"]

        if pdb_ids:
            selected_pdb = st.selectbox(
                "Select a PDB structure",
                pdb_ids,
                key="selected_pdb"
            )
            st.markdown(
                f"[Open interactive 3D structure]"
                f"(https://www.rcsb.org/3d-view/{selected_pdb})"
            )

            pdb_file = download_pdb(selected_pdb)
            if pdb_file:
                st.download_button(
                    "Download PDB file",
                    data=pdb_file,
                    file_name=f"{selected_pdb}.pdb",
                    mime="chemical/x-pdb"
                )
            else:
                st.warning("The PDB file could not be downloaded.")
        else:
            st.warning(
                "No structures linked to this UniProt accession were returned."
            )

    st.caption(
        "A PDB hit must be checked for species, chain, construct and sequence "
        "coverage before claiming that it represents the complete human protein."
    )

    st.divider()
    st.caption(
        "Bioinformatics Dashboard | Gene & Protein Characterization"
    )
