
"""Gene-independent transcript and isoform analysis for Helixora Bioinformatics Studio.

This module retrieves transcript records through NCBI Gene -> NCBI Nucleotide,
annotates GenBank features, compares transcript sequences, reports annotated
exon/CDS/protein information, and transfers selected records to Primer Design.

Important interpretation limits:
- GenBank feature annotations are reported as annotated, not experimentally proven.
- Exact-match distinguishing windows are compared only against the retrieved set.
- This module does not design primers or establish genome-wide specificity.
- Alternative-splicing event labels are not inferred unless exon/junction evidence
  is explicitly available in the source annotations.
"""

from __future__ import annotations

import re
from collections import Counter
from typing import Iterable

import pandas as pd
import streamlit as st

from Bio import Entrez, SeqIO
from Bio.Seq import Seq

from core.data_manager import register_records


# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

DEFAULT_GENE = "MEF2C"
DEFAULT_ORGANISM = "Homo sapiens"
DEFAULT_LIMIT = 250
DEFAULT_WINDOW_SIZE = 22
BATCH_SIZE = 50


# -----------------------------------------------------------------------------
# NCBI GENE RESOLUTION AND TRANSCRIPT RETRIEVAL
# -----------------------------------------------------------------------------

def _validate_search_input(
    gene_query: str,
    organism: str
) -> tuple[str, str]:

    gene_query = str(gene_query or "").strip()
    organism = str(organism or "").strip()

    if not gene_query:
        raise ValueError("Enter a gene symbol or NCBI Gene ID.")

    if not organism:
        raise ValueError("Enter an organism name.")

    return gene_query, organism


def _extract_gene_summaries(summaries):
    """Normalize different Biopython Entrez Gene ESummary response structures."""

    if isinstance(summaries, dict):
        document_set = summaries.get(
            "DocumentSummarySet",
            summaries
        )

        if isinstance(document_set, dict):
            summary_items = document_set.get(
                "DocumentSummary",
                []
            )
        else:
            summary_items = document_set

    else:
        summary_items = summaries

    if summary_items is None:
        return []

    if not isinstance(summary_items, (list, tuple)):
        summary_items = [summary_items]

    return [
        item for item in summary_items
        if hasattr(item, "get")
    ]


def resolve_gene_id(
    gene_query: str,
    organism: str
) -> tuple[str, str]:
    """Resolve a gene symbol to an NCBI Gene ID."""

    gene_query, organism = _validate_search_input(
        gene_query,
        organism
    )

    is_numeric_id = bool(
        re.fullmatch(r"\d+", gene_query)
    )

    if is_numeric_id:
        term = (
            f"({gene_query}[GeneID]) "
            f'AND ("{organism}"[Organism])'
        )
    else:
        safe_symbol = gene_query.replace('"', "").strip()
        term = (
            f'("{safe_symbol}"[Gene Name]) '
            f'AND ("{organism}"[Organism])'
        )

    # Search NCBI Gene
    with Entrez.esearch(
        db="gene",
        term=term,
        retmax=20,
        sort="relevance"
    ) as handle:
        search_result = Entrez.read(handle)

    gene_ids = [
        str(gid)
        for gid in search_result.get("IdList", [])
    ]

    if not gene_ids:
        raise LookupError(
            f"No NCBI Gene record was found for "
            f"'{gene_query}' in {organism}."
        )

    # Numeric Gene ID validation
    if is_numeric_id:
        if gene_query not in gene_ids:
            raise LookupError(
                "The supplied Gene ID does not match "
                "the selected organism."
            )
        return gene_query, gene_query

    # A single result from the exact Gene Name + Organism query
    # is sufficient to resolve the gene.
    if len(gene_ids) == 1:
        return gene_ids[0], gene_query

    # Multiple results: inspect Gene summaries
    with Entrez.esummary(
        db="gene",
        id=",".join(gene_ids),
        retmode="xml"
    ) as handle:
        summaries = Entrez.read(handle)

    summary_items = _extract_gene_summaries(summaries)

    wanted = gene_query.casefold()
    exact = []

    for item in summary_items:

        symbol = str(
            item.get("Name", "")
        ).strip()

        gene_id = str(
            item.get("Id")
            or item.get("uid")
            or item.get("GeneId")
            or ""
        ).strip()

        if symbol.casefold() == wanted and gene_id:
            exact.append((gene_id, symbol))

    if len(exact) == 1:
        return exact[0]

    if len(exact) > 1:
        choices = ", ".join(
            gid for gid, _ in exact
        )
        raise LookupError(
            f"Multiple exact matches found for '{gene_query}'. "
            f"Please use a Gene ID: {choices}."
        )

    raise LookupError(
        f"NCBI returned multiple candidate records for "
        f"'{gene_query}', but an exact match could not be confirmed. "
        f"Candidate Gene IDs: {', '.join(gene_ids[:10])}. "
        "Please enter the intended NCBI Gene ID."
    )

def build_nuccore_query(
    gene_id: str,
    organism: str
) -> str:
    """Build a nucleotide query using a resolved NCBI Gene ID."""

    gene_id = str(gene_id).strip()
    organism = str(organism).strip()

    if not gene_id.isdigit():
        raise ValueError(
            "A numeric NCBI Gene ID is required."
        )

    if not organism:
        raise ValueError(
            "An organism name is required."
        )

    return (
        f"({gene_id}[GeneID]) "
        f'AND ("{organism}"[Organism]) '
        f"AND (biomol_mrna[PROP])"
    )


def build_gene_query(
    gene_query: str,
    organism: str
) -> str:
    """Backward-compatible query builder."""

    gene_id, _ = resolve_gene_id(
        gene_query,
        organism
    )

    return build_nuccore_query(
        gene_id,
        organism
    )


def fetch_gene_transcripts(
    gene_query: str,
    organism: str,
    limit: int = DEFAULT_LIMIT
) -> tuple[list, int, str]:
    """Resolve a gene, retrieve GenBank mRNA records, and return records/count/GeneID."""

    gene_query, organism = _validate_search_input(
        gene_query,
        organism
    )

    limit = max(
        1,
        min(int(limit), 5000)
    )

    gene_id, _resolved_name = resolve_gene_id(
        gene_query,
        organism
    )

    term = build_nuccore_query(
        gene_id,
        organism
    )

    # Search NCBI Nucleotide
    with Entrez.esearch(
        db="nuccore",
        term=term,
        retmax=limit,
        sort="relevance",
        usehistory="n"
    ) as handle:

        search_result = Entrez.read(handle)

    ids = list(
        search_result.get("IdList", [])
    )

    total = int(
        search_result.get("Count", 0)
    )

    if not ids:
        return [], total, gene_id

    # Fetch records in batches
    records = []

    for offset in range(
        0,
        len(ids),
        BATCH_SIZE
    ):

        batch_ids = ids[
            offset:offset + BATCH_SIZE
        ]

        with Entrez.efetch(
            db="nuccore",
            id=",".join(batch_ids),
            rettype="gb",
            retmode="text"
        ) as handle:

            batch_records = list(
                SeqIO.parse(handle, "genbank")
            )

        records.extend(batch_records)

    # Deduplicate by accession while preserving NCBI relevance order
    unique = {}

    for record in records:

        unique.setdefault(
            record.id,
            record
        )

    return (
        list(unique.values()),
        total,
        gene_id
    )


# -----------------------------------------------------------------------------
# GENBANK ANNOTATION HELPERS
# -----------------------------------------------------------------------------

def get_qualifier(
    feature,
    key: str
) -> str:

    if feature is None:
        return ""

    values = feature.qualifiers.get(
        key,
        []
    )

    return "; ".join(
        str(value) for value in values
    ) if values else ""


def get_cds_feature(record):

    return next(
        (
            feature
            for feature in record.features
            if feature.type == "CDS"
        ),
        None
    )


def get_exon_features(record) -> list:

    return [
        feature
        for feature in record.features
        if feature.type == "exon"
    ]


def get_exon_count(record):

    features = get_exon_features(record)

    return len(features) if features else None


def classify_accession(
    accession: str
) -> str:

    accession = str(accession or "")

    if accession.startswith("NM_"):
        return "RefSeq curated mRNA"

    if accession.startswith("NR_"):
        return "RefSeq curated non-coding RNA"

    if accession.startswith("XM_"):
        return "RefSeq predicted mRNA"

    if accession.startswith("XR_"):
        return "RefSeq predicted non-coding RNA"

    if accession.startswith("ENST"):
        return "Ensembl transcript"

    return "Other nucleotide record"


def _protein_details(record) -> tuple[str, str]:

    protein_ids = set()
    protein_lengths = set()

    for feature in record.features:

        if feature.type != "CDS":
            continue

        protein_id = get_qualifier(
            feature,
            "protein_id"
        )

        translation = get_qualifier(
            feature,
            "translation"
        )

        if protein_id:

            protein_ids.update(
                x.strip()
                for x in protein_id.split(";")
                if x.strip()
            )

        if translation:

            protein_lengths.add(
                str(
                    len(
                        translation.replace("*", "").replace(" ", "")
                    )
                )
            )

    return (
        "; ".join(sorted(protein_ids)) or "—",
        "; ".join(
            sorted(
                protein_lengths,
                key=lambda x: int(x)
            )
        ) or "—"
    )


def _feature_location(feature) -> str:

    try:
        return str(feature.location)

    except Exception:
        return "Not available"


def catalogue(records) -> pd.DataFrame:
    """Create a transcript catalogue using fields supported by GenBank records."""

    rows = []

    for record in records or []:

        cds = get_cds_feature(record)

        protein_accession, protein_length = (
            _protein_details(record)
        )

        gene_features = [
            feature
            for feature in record.features
            if feature.type == "gene"
        ]

        gene_names = sorted({
            get_qualifier(feature, "gene")
            for feature in gene_features
            if get_qualifier(feature, "gene")
        })

        exon_features = get_exon_features(record)

        cds_length = (
            len(cds.extract(record.seq))
            if cds is not None
            else 0
        )

        rows.append({

            "Accession": record.id,

            "Record class": classify_accession(
                record.id
            ),

            "Description": record.description,

            "Organism": record.annotations.get(
                "organism",
                "Not annotated"
            ),

            "Transcript length (nt)": len(record.seq),

            # Keep this column numeric; use None for missing annotations
            # so Streamlit/PyArrow does not receive mixed int/string values.
            "Exon count": (
                len(exon_features)
                if exon_features
                else None
            ),

            "CDS length (nt)": (
                cds_length
                if cds is not None
                else None
            ),

            "CDS location": (
                _feature_location(cds)
                if cds is not None
                else "Not annotated"
            ),

            "Protein accession": protein_accession,

            "Protein length (aa)": protein_length,

            "Coding status": (
                "CDS annotated"
                if cds is not None
                else "No CDS feature"
            ),

            "Gene annotation": (
                "; ".join(gene_names)
                if gene_names
                else "Not annotated"
            ),

            "Annotation evidence": "GenBank feature annotation"

        })

    return pd.DataFrame(rows)


def normalize_sequence(sequence) -> str:

    return (
        re.sub(
            r"\s+",
            "",
            str(sequence).upper()
        )
        .replace("U", "T")
    )


def gc_content(sequence: str) -> float:

    seq = normalize_sequence(sequence)

    valid = [
        base
        for base in seq
        if base in "ACGT"
    ]

    if not valid:
        return 0.0

    return round(
        100 * (
            valid.count("G") + valid.count("C")
        ) / len(valid),
        2
    )


def compare_transcripts(records) -> pd.DataFrame:

    rows = []

    for record in records or []:

        seq = normalize_sequence(record.seq)

        cds = get_cds_feature(record)

        exon_features = get_exon_features(record)

        rows.append({

            "Accession": record.id,

            "Transcript length (nt)": len(seq),

            "GC content (%)": gc_content(seq),

            # Keep this column numeric; use None for missing annotations
            # so Streamlit/PyArrow does not receive mixed int/string values.
            "Exon count": (
                len(exon_features)
                if exon_features
                else None
            ),

            "CDS length (nt)": (
                len(cds.extract(record.seq))
                if cds is not None
                else None
            ),

            "CDS annotated": (
                "Yes"
                if cds is not None
                else "No"
            ),

            "Strand": (
                str(cds.location.strand)
                if cds is not None
                else "Not annotated"
            )

        })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# TRANSCRIPT-DISCRIMINATING REGIONS
# NOT PRIMER DESIGN
# -----------------------------------------------------------------------------

def find_unique_windows(
    target: str,
    other_sequences: Iterable[str],
    k: int = DEFAULT_WINDOW_SIZE
):

    target = normalize_sequence(target)

    others = [
        normalize_sequence(seq)
        for seq in other_sequences
    ]

    if k < 1 or len(target) < k:
        return []

    windows = []
    seen = set()

    for start in range(
        len(target) - k + 1
    ):

        window = target[
            start:start + k
        ]

        if not re.fullmatch(
            r"[ACGT]+",
            window
        ):
            continue

        if window in seen:
            continue

        reverse_complement = str(
            Seq(window).reverse_complement()
        )

        if all(
            window not in seq
            and reverse_complement not in seq
            for seq in others
        ):

            windows.append({

                "Start (1-based)": start + 1,

                "End (1-based)": start + k,

                "Sequence (5′→3′)": window,

                "Length (nt)": k

            })

            seen.add(window)

    return windows


def transcript_discriminating_regions(
    records,
    k: int = DEFAULT_WINDOW_SIZE
) -> pd.DataFrame:

    sequences = [
        normalize_sequence(record.seq)
        for record in records
    ]

    rows = []

    for index, record in enumerate(records):

        others = (
            sequences[:index]
            + sequences[index + 1:]
        )

        windows = find_unique_windows(
            sequences[index],
            others,
            k=k
        )

        rows.append({

            "Accession": record.id,

            "Unique exact-match windows": len(windows),

            "Example region": (
                windows[0]["Sequence (5′→3′)"]
                if windows
                else "—"
            ),

            "Interpretation": (

                "Distinguishing sequence window(s) found "
                "within the retrieved set"

                if windows

                else

                "No distinguishing exact-match window "
                "found within the retrieved set"

            ),

            "Evidence level": (
                "Computationally inferred; retrieved records only"
            )

        })

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# FUNCTIONAL INTERPRETATION AND EXPORT
# -----------------------------------------------------------------------------

def render_functional_evidence_section():

    st.subheader(
        "Functional Isoform Interpretation"
    )

    st.markdown(
        "Functional differences between isoforms require accession-specific "
        "evidence. Transcript variant numbering alone does not establish "
        "whether a protein is an activator, repressor, or has a distinct "
        "biological function. Consider curated protein annotations, domain "
        "composition, exon inclusion, interaction studies, and gene-specific "
        "functional experiments."
    )

    st.info(
        "This module does not automatically assign functional labels "
        "without supporting evidence."
    )


def _download_fasta(
    records,
    filename: str,
    label: str
):

    fasta = "".join(
        record.format("fasta")
        for record in records
    )

    st.download_button(
        label,
        data=fasta,
        file_name=filename,
        mime="text/plain"
    )


def _reset_transcript_results():

    for key in (
        "isoform_records",
        "isoform_total",
        "isoform_gene_id",
        "isoform_feasibility"
    ):

        st.session_state.pop(
            key,
            None
        )


# -----------------------------------------------------------------------------
# STREAMLIT PAGE
# -----------------------------------------------------------------------------

def render_transcript_analysis_page():

    st.title(
        "Transcript & Isoform Analysis"
    )

    st.caption(
        "Gene identification • Transcript annotation • "
        "Sequence comparison • Isoform evidence"
    )

    st.write(
        "Retrieve transcript records for a gene, inspect their annotations, "
        "compare sequence features, and select records for downstream analysis."
    )

    st.info(
        "NCBI requests use the shared Entrez configuration "
        "provided by the application."
    )

    # -------------------------------------------------------------------------
    # 1. Gene identification and transcript retrieval
    # -------------------------------------------------------------------------

    st.subheader(
        "1. Gene Identification and Transcript Retrieval"
    )

    col1, col2, col3 = st.columns(
        [1.2, 1.5, 1]
    )

    with col1:

        gene_query = st.text_input(
            "Gene symbol or NCBI Gene ID",
            value=st.session_state.get(
                "isoform_gene_query",
                DEFAULT_GENE
            ),
            placeholder="e.g. TP53, BRCA1, 7157",
            key="transcript_gene_input"
        )

    with col2:

        organism = st.text_input(
            "Organism",
            value=st.session_state.get(
                "isoform_organism",
                DEFAULT_ORGANISM
            ),
            placeholder="Homo sapiens",
            key="transcript_organism_input"
        )

    with col3:

        limit = st.number_input(
            "Maximum transcript records",
            min_value=1,
            max_value=1000,
            value=DEFAULT_LIMIT,
            step=10
        )

    if st.button(
        "Retrieve Transcripts from NCBI",
        type="primary",
        use_container_width=True
    ):

        try:

            with st.spinner(
                "Resolving gene through NCBI Gene "
                "and retrieving transcript records..."
            ):

                records, total, gene_id = fetch_gene_transcripts(
                    gene_query,
                    organism,
                    int(limit)
                )

            _reset_transcript_results()

            st.session_state["isoform_records"] = records

            st.session_state["isoform_total"] = total

            st.session_state["isoform_gene_id"] = gene_id

            st.session_state["isoform_gene_query"] = (
                gene_query.strip()
            )

            st.session_state["isoform_organism"] = (
                organism.strip()
            )

            if records:

                register_records(
                    records,
                    source="Transcript & Isoform Analysis",
                    gene=gene_query.strip(),
                    organism=organism.strip()
                )

                st.success(
                    f"Retrieved {len(records)} records. "
                    f"NCBI reported {total} matching nucleotide records."
                )

            else:

                st.warning(
                    "The gene was resolved, but no mRNA records "
                    "were returned for this organism."
                )

        except Exception as exc:

            st.error(
                f"NCBI retrieval failed: {exc}"
            )

    # -------------------------------------------------------------------------
    # Load retrieved records
    # -------------------------------------------------------------------------

    records = st.session_state.get(
        "isoform_records",
        []
    )

    if not records:

        st.info(
            "Enter a gene and organism, then retrieve transcript records."
        )

        return

    active_gene = st.session_state.get(
        "isoform_gene_query",
        gene_query
    )

    safe_gene = re.sub(
        r"[^A-Za-z0-9_-]+",
        "_",
        active_gene
    ).strip("_") or "gene"

    df = catalogue(records)

    # -------------------------------------------------------------------------
    # 2. Transcript catalogue
    # -------------------------------------------------------------------------

    st.divider()

    st.subheader(
        "2. Transcript Catalogue"
    )

    m1, m2, m3, m4 = st.columns(4)

    m1.metric(
        "Retrieved records",
        len(records)
    )

    m2.metric(
        "Curated RefSeq",
        int(
            df["Record class"].str.contains(
                "curated",
                case=False
            ).sum()
        )
    )

    m3.metric(
        "Predicted RefSeq",
        int(
            df["Record class"].str.contains(
                "predicted",
                case=False
            ).sum()
        )
    )

    m4.metric(
        "CDS annotated",
        int(
            df["Coding status"].eq(
                "CDS annotated"
            ).sum()
        )
    )

    st.caption(
        f"Resolved NCBI Gene ID: "
        f"{st.session_state.get('isoform_gene_id', 'Not available')} "
        f"| NCBI nucleotide count: "
        f"{st.session_state.get('isoform_total', 'Not available')}"
    )

    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True
    )

    st.download_button(
        "Download Transcript Catalogue (CSV)",
        df.to_csv(index=False).encode("utf-8-sig"),
        file_name=f"{safe_gene}_transcript_catalogue.csv",
        mime="text/csv"
    )

    _download_fasta(
        records,
        f"{safe_gene}_transcripts.fasta",
        "Download Transcript Sequences (FASTA)"
    )

    # -------------------------------------------------------------------------
    # 3. Transcript-to-protein mapping
    # -------------------------------------------------------------------------

    st.divider()

    st.subheader(
        "3. Transcript-to-Protein Mapping"
    )

    protein_df = df.loc[
        df["Protein accession"].ne("—"),
        [
            "Accession",
            "Protein accession",
            "Protein length (aa)",
            "CDS length (nt)"
        ]
    ]

    if protein_df.empty:

        st.warning(
            "No annotated protein accession was found "
            "in the retrieved records."
        )

    else:

        st.dataframe(
            protein_df,
            use_container_width=True,
            hide_index=True
        )

    # -------------------------------------------------------------------------
    # 4. Transcript structural and sequence comparison
    # -------------------------------------------------------------------------

    st.divider()

    st.subheader(
        "4. Transcript Structural and Sequence Comparison"
    )

    comparison = compare_transcripts(
        records
    )

    st.dataframe(
        comparison,
        use_container_width=True,
        hide_index=True
    )

    st.download_button(
        "Download Transcript Comparison (CSV)",
        comparison.to_csv(index=False).encode("utf-8-sig"),
        file_name=f"{safe_gene}_transcript_comparison.csv",
        mime="text/csv"
    )

    # -------------------------------------------------------------------------
    # 5. Transcript-discriminating region analysis
    # -------------------------------------------------------------------------

    st.divider()

    st.subheader(
        "5. Transcript-Discriminating Region Analysis"
    )

    st.caption(
        "This is an exact-match comparison among retrieved transcript "
        "sequences, not primer design or genome-wide specificity testing."
    )

    window_size = st.select_slider(
        "Sequence window length (nt)",
        options=[18, 20, 22, 24, 26, 28],
        value=DEFAULT_WINDOW_SIZE
    )

    if st.button(
        "Analyze Transcript-Discriminating Regions"
    ):

        try:

            st.session_state["isoform_feasibility"] = (
                transcript_discriminating_regions(
                    records,
                    int(window_size)
                )
            )

        except Exception as exc:

            st.error(
                f"Region analysis failed: {exc}"
            )

    region_df = st.session_state.get(
        "isoform_feasibility"
    )

    if region_df is not None:

        st.dataframe(
            region_df,
            use_container_width=True,
            hide_index=True
        )

        st.download_button(
            "Download Region Analysis (CSV)",
            region_df.to_csv(index=False).encode("utf-8-sig"),
            file_name=f"{safe_gene}_discriminating_regions.csv",
            mime="text/csv"
        )

    # -------------------------------------------------------------------------
    # 6. Functional isoform evidence
    # -------------------------------------------------------------------------

    st.divider()

    st.subheader(
        "6. Functional Isoform Evidence"
    )

    render_functional_evidence_section()

    # -------------------------------------------------------------------------
    # 7. Transfer selected transcripts to Primer Design
    # -------------------------------------------------------------------------

    st.divider()

    st.subheader(
        "7. Transfer Selected Transcripts to Primer Design"
    )

    accessions = df[
        "Accession"
    ].tolist()

    curated = df.loc[
        df["Record class"].str.contains(
            "curated",
            case=False
        ),
        "Accession"
    ].tolist()

    default_selection = (
        curated or accessions
    )[:10]

    selected = st.multiselect(
        "Select transcript records",
        options=accessions,
        default=default_selection,
        key="transcript_records_to_transfer"
    )

    if st.button(
        "Send Selected Transcripts to Primer Design",
        type="primary"
    ):

        selected_set = set(selected)

        selected_records = [
            record
            for record in records
            if record.id in selected_set
        ]

        if not selected_records:

            st.warning(
                "Select at least one transcript record."
            )

        else:

            st.session_state["primer_transcripts"] = (
                selected_records
            )

            st.session_state["primer_search_gene"] = (
                active_gene
            )

            st.session_state["primer_search_organism"] = (
                st.session_state.get(
                    "isoform_organism",
                    organism
                )
            )

            # Invalidate any alignment tied to a previous selection.
            for key in (
                "primer_alignment_text",
                "primer_alignment_fasta",
                "primer_alignment_job_id",
                "primer_alignment_summary",
                "primer_conservation_df"
            ):

                st.session_state[key] = None

            register_records(
                selected_records,
                source="Transcript Analysis selection",
                gene=active_gene,
                organism=st.session_state.get(
                    "isoform_organism",
                    organism
                )
            )

            st.success(
                f"{len(selected_records)} transcript(s) "
                "transferred to Primer Design."
            )

    # -------------------------------------------------------------------------
    # 8. Experimental validation considerations
    # -------------------------------------------------------------------------

    st.divider()

    st.subheader(
        "8. Experimental Validation Considerations"
    )

    st.markdown(
        "- Review transcript annotations and exon structures before selecting an assay target.\n"
        "- Confirm transcript-discriminating regions against the current reference genome.\n"
        "- Design primer pairs in the Primer Design module.\n"
        "- Use an appropriate genome-wide specificity workflow "
        "(for example, NCBI Primer-BLAST) for final in-silico validation.\n"
        "- Include no-template and no-RT controls where appropriate.\n"
        "- Confirm product identity and amplification efficiency experimentally."
    )


# -----------------------------------------------------------------------------
# BACKWARD COMPATIBILITY
# -----------------------------------------------------------------------------

def render_mef2c_page():

    render_transcript_analysis_page()