"""
Primer Design Module
====================

This module provides three input routes for primer design:

1. Shared transcript registry
2. Direct NCBI retrieval
3. FASTA file upload

Primer generation and evaluation are delegated to:
modules.rtqpcr_primer

The module is designed to work with:
core.data_manager
modules.transcript_analysis
"""

from __future__ import annotations

import io
from typing import Any, Iterable

import streamlit as st
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from modules.ncbi_client import (
    search_transcripts,
    fetch_genbank_record,
)

from modules.rtqpcr_primer import (
    render_rtqpcr_primer_design,
)

from core.data_manager import register_records

try:
    from core.data_manager import get_registered_records
except ImportError:
    get_registered_records = None


# ============================================================
# ORGANISM DATABASE
# ============================================================

COMMON_ORGANISMS = [
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Danio rerio",
    "Drosophila melanogaster",
    "Caenorhabditis elegans",
    "Gallus gallus",
    "Bos taurus",
    "Sus scrofa",
    "Ovis aries",
    "Canis lupus familiaris",
    "Macaca mulatta",
    "Pan troglodytes",
    "Xenopus laevis",
    "Xenopus tropicalis",
    "Oryctolagus cuniculus",
    "Mesocricetus auratus",
    "Cavia porcellus",
    "Arabidopsis thaliana",
    "Zea mays",
    "Oryza sativa",
    "Triticum aestivum",
    "Saccharomyces cerevisiae",
    "Escherichia coli",
    "Bacillus subtilis",
    "Other (Enter manually)",
]


# ============================================================
# SESSION STATE
# ============================================================

def initialize_primer_state():

    defaults = {

        "primer_candidates": [],

        "primer_transcripts": [],

        "primer_search_gene": "",

        "primer_search_organism": "Homo sapiens",

        "primer_source_mode": "Shared transcript registry",

        "primer_registry_signature": None,

    }

    for key, value in defaults.items():

        if key not in st.session_state:

            st.session_state[key] = value


# ============================================================
# RECORD UTILITIES
# ============================================================

def get_record_id(record):

    return str(
        getattr(record, "id", "")
        or getattr(record, "name", "")
        or "unknown"
    )


def get_record_description(record):

    return str(
        getattr(record, "description", "")
        or get_record_id(record)
    )


def get_record_sequence(record):

    return "".join(
        str(getattr(record, "seq", "")).split()
    ).upper()


def validate_nucleotide_sequence(sequence):

    if not sequence:

        return False

    allowed = set(
        "ACGTUNRYKMSWBDHV"
    )

    return not (
        set(sequence) - allowed
    )


def normalize_records(records):

    """
    Normalize sequence records and remove duplicates.

    Only valid nucleotide sequences are accepted.
    """

    normalized = []

    seen = set()

    for item in records or []:

        try:

            sequence = get_record_sequence(item)

            record_id = get_record_id(item).strip()

            if not record_id:

                continue

            if not validate_nucleotide_sequence(sequence):

                continue

            key = record_id.casefold()

            if key in seen:

                continue

            seen.add(key)

            record = SeqRecord(

                Seq(sequence),

                id=record_id,

                name=str(
                    getattr(item, "name", record_id)
                    or record_id
                ),

                description=get_record_description(item),

                annotations=dict(
                    getattr(item, "annotations", {}) or {}
                ),

                dbxrefs=list(
                    getattr(item, "dbxrefs", []) or []
                ),

            )

            normalized.append(record)

        except Exception:

            continue

    return normalized


# ============================================================
# SHARED REGISTRY
# ============================================================

def get_shared_registry_records():

    """
    Retrieve records from the shared transcript registry.
    """

    if get_registered_records is None:

        return []

    try:

        registered = get_registered_records()

        if isinstance(registered, dict):

            registered = list(
                registered.values()
            )

        records = []

        for item in registered or []:

            if isinstance(item, dict):

                record = item.get("record")

            else:

                record = item

            if record is not None:

                records.append(record)

        return normalize_records(records)

    except Exception:

        return []


def create_registry_signature(records):

    return tuple(

        (
            get_record_id(record),
            len(record.seq)
        )

        for record in records

    )


def synchronize_registry():

    """
    Synchronize shared registry records with local primer state.
    """

    records = get_shared_registry_records()

    signature = create_registry_signature(records)

    previous_signature = (
        st.session_state.primer_registry_signature
    )

    if records and signature != previous_signature:

        st.session_state.primer_transcripts = records

        st.session_state.primer_registry_signature = signature


def register_sequence_records(
    records,
    gene="",
    organism="",
    source="Primer Design",
):

    if not records:

        return

    try:

        register_records(

            records,

            source=source,

            gene=gene,

            organism=organism,

        )

    except TypeError:

        # Compatibility with older data_manager versions.

        register_records(
            records,
            source=source
        )


# ============================================================
# NCBI SEARCH
# ============================================================

def search_gene_transcripts(
    gene,
    organism,
    limit,
):

    return search_transcripts(

        gene_symbol=gene.strip(),

        organism=organism.strip(),

        limit=int(limit),

    )


# ============================================================
# NCBI RETRIEVAL
# ============================================================

def retrieve_selected_transcripts(
    selected_candidates
):

    records = []

    errors = []

    for candidate in selected_candidates:

        accession = (

            candidate.get("accession")

            or candidate.get("id")

            or candidate.get("uid")

        )

        if not accession:

            errors.append(
                "Candidate without accession."
            )

            continue

        try:

            record = fetch_genbank_record(
                str(accession)
            )

            records.append(record)

        except Exception as error:

            errors.append(
                f"{accession}: {error}"
            )

    records = normalize_records(records)

    if errors:

        st.warning(
            "Some records could not be retrieved:\n\n"
            + "\n".join(errors)
        )

    return records


# ============================================================
# FASTA EXPORT
# ============================================================

def export_fasta(records):

    buffer = io.StringIO()

    valid_records = normalize_records(records)

    SeqIO.write(
        valid_records,
        buffer,
        "fasta"
    )

    return buffer.getvalue()


# ============================================================
# ORGANISM SELECTOR
# ============================================================

def organism_selector():

    current_organism = (
        st.session_state.primer_search_organism
    )

    if current_organism in COMMON_ORGANISMS:

        default_index = COMMON_ORGANISMS.index(
            current_organism
        )

    else:

        default_index = (
            len(COMMON_ORGANISMS) - 1
        )

    selected = st.selectbox(

        "Organism",

        options=COMMON_ORGANISMS,

        index=default_index,

        key="primer_organism_select",

    )

    if selected == "Other (Enter manually)":

        organism = st.text_input(

            "Scientific name",

            value=(
                current_organism
                if current_organism not in COMMON_ORGANISMS
                else ""
            ),

            placeholder="Example: Rattus rattus",

            key="primer_custom_organism",

        )

        return organism.strip()

    return selected


# ============================================================
# CANDIDATE LABEL
# ============================================================

def candidate_label(candidate):

    return str(

        candidate.get("label")

        or candidate.get("accession")

        or candidate.get("id")

        or candidate.get("uid")

        or "Unknown transcript"

    )


# ============================================================
# SHARED REGISTRY INTERFACE
# ============================================================

def render_shared_registry():

    st.subheader(
        "Transcript Analysis Integration"
    )

    records = get_shared_registry_records()

    if not records:

        st.info(
            """
            No transcript records are currently available.

            Recommended workflow:

            1. Open Transcript & Isoform Analysis.
            2. Retrieve the desired gene transcripts.
            3. Select the transcript records.
            4. Transfer them to Primer Design.

            You can also retrieve sequences directly from NCBI
            or upload a FASTA file.
            """
        )

        return

    st.success(
        f"{len(records)} transcript record(s) found."
    )

    table = pd.DataFrame([

        {

            "Accession": get_record_id(record),

            "Description": get_record_description(record),

            "Length (nt)": len(record.seq),

        }

        for record in records

    ])

    st.dataframe(

        table,

        use_container_width=True,

        hide_index=True,

    )

    selected_ids = st.multiselect(

        "Select transcripts for primer design",

        options=[
            get_record_id(record)
            for record in records
        ],

        default=[
            get_record_id(record)
            for record in records
        ],

        key="primer_registry_selection",

    )

    if st.button(

        "Use selected registry transcripts",

        key="use_registry_records",

    ):

        selected_records = [

            record

            for record in records

            if get_record_id(record) in selected_ids

        ]

        if not selected_records:

            st.warning(
                "Select at least one transcript."
            )

            return

        st.session_state.primer_transcripts = (
            selected_records
        )

        st.session_state.primer_registry_signature = (
            create_registry_signature(records)
        )

        st.session_state.primer_source_mode = (
            "Shared transcript registry"
        )

        st.success(
            f"{len(selected_records)} transcript(s) selected."
        )


# ============================================================
# DIRECT NCBI INTERFACE
# ============================================================

def render_ncbi_search():

    st.subheader(
        "Direct NCBI Transcript Retrieval"
    )

    col1, col2 = st.columns(2)

    with col1:

        gene = st.text_input(

            "Gene symbol",

            value=st.session_state.primer_search_gene,

            key="primer_gene_input",

            placeholder="Example: TP53",

        )

    with col2:

        organism = organism_selector()

    limit = st.number_input(

        "Maximum transcript candidates",

        min_value=1,

        max_value=100,

        value=20,

        step=1,

        key="primer_candidate_limit",

    )

    if st.button(

        "Search NCBI",

        type="primary",

        key="primer_ncbi_search",

    ):

        if not gene.strip() or not organism.strip():

            st.error(
                "Gene symbol and organism are required."
            )

            return

        try:

            with st.spinner(
                "Searching NCBI..."
            ):

                candidates = search_gene_transcripts(

                    gene,

                    organism,

                    int(limit),

                )

            st.session_state.primer_candidates = (
                candidates or []
            )

            st.session_state.primer_search_gene = (
                gene.strip()
            )

            st.session_state.primer_search_organism = (
                organism.strip()
            )

            if candidates:

                st.success(
                    f"{len(candidates)} candidates found."
                )

            else:

                st.warning(
                    "No transcript candidates found."
                )

        except Exception as error:

            st.error(
                f"NCBI search failed: {error}"
            )

    candidates = st.session_state.primer_candidates

    if not candidates:

        return

    labels = [
        candidate_label(candidate)
        for candidate in candidates
    ]

    selected_labels = st.multiselect(

        "Select transcript candidates",

        options=labels,

        key="primer_ncbi_candidates_selected",

    )

    if st.button(

        "Download selected GenBank records",

        key="primer_download_ncbi",

    ):

        selected_candidates = [

            candidate

            for candidate in candidates

            if candidate_label(candidate) in selected_labels

        ]

        if not selected_candidates:

            st.warning(
                "Select at least one candidate."
            )

            return

        with st.spinner(
            "Downloading transcript records..."
        ):

            records = retrieve_selected_transcripts(
                selected_candidates
            )

        if records:

            st.session_state.primer_transcripts = records

            st.session_state.primer_source_mode = "NCBI"

            register_sequence_records(

                records,

                gene=st.session_state.primer_search_gene,

                organism=st.session_state.primer_search_organism,

                source="Primer Design / NCBI",

            )

            st.success(
                f"{len(records)} records downloaded."
            )


# ============================================================
# FASTA UPLOAD INTERFACE
# ============================================================

def render_fasta_upload():

    st.subheader(
        "FASTA Sequence Upload"
    )

    uploaded_file = st.file_uploader(

        "Upload nucleotide FASTA file",

        type=[
            "fa",
            "fasta",
            "fna",
            "txt",
        ],

        key="primer_fasta_upload",

    )

    if uploaded_file is None:

        return

    try:

        content = uploaded_file.getvalue().decode(
            "utf-8-sig"
        )

        parsed_records = list(

            SeqIO.parse(

                io.StringIO(content),

                "fasta"

            )

        )

        records = normalize_records(
            parsed_records
        )

        if not records:

            st.error(
                "No valid nucleotide sequences were found."
            )

            return

        st.success(
            f"{len(records)} valid sequences detected."
        )

        table = pd.DataFrame([

            {

                "Sequence ID": record.id,

                "Description": record.description,

                "Length (nt)": len(record.seq),

            }

            for record in records

        ])

        st.dataframe(

            table,

            use_container_width=True,

            hide_index=True,

        )

        if st.button(

            "Use uploaded FASTA sequences",

            key="primer_use_uploaded_fasta",

        ):

            st.session_state.primer_transcripts = records

            st.session_state.primer_source_mode = (
                "FASTA upload"
            )

            register_sequence_records(

                records,

                source="Primer Design / FASTA",

            )

            st.success(
                "FASTA sequences are ready for primer design."
            )

    except Exception as error:

        st.error(
            f"FASTA parsing failed: {error}"
        )


# ============================================================
# ACTIVE TRANSCRIPT TABLE
# ============================================================

def render_active_transcripts(records):

    st.subheader(
        "Active Primer Design Targets"
    )

    rows = []

    for record in records:

        sequence = str(record.seq).upper()

        length = len(sequence)

        gc_count = (
            sequence.count("G")
            + sequence.count("C")
        )

        gc_percentage = (

            100 * gc_count / length

            if length > 0

            else 0

        )

        rows.append({

            "Accession / ID": record.id,

            "Description": record.description,

            "Length (nt)": length,

            "GC (%)": round(
                gc_percentage,
                2
            ),

        })

    table = pd.DataFrame(rows)

    st.dataframe(

        table,

        use_container_width=True,

        hide_index=True,

    )

    st.download_button(

        "Download active sequences (FASTA)",

        data=export_fasta(records),

        file_name="primer_design_input_transcripts.fasta",

        mime="text/plain",

        key="primer_download_active_fasta",

    )


# ============================================================
# MAIN PAGE
# ============================================================

def render_primer_design_page():

    initialize_primer_state()

    synchronize_registry()

    st.title(
        "🧬 Primer Design"
    )

    st.markdown(
        """
        Design and evaluate PCR and RT-qPCR primers
        using verified nucleotide transcript sequences.
        """
    )

    st.info(
        """
        Recommended workflow:

        Retrieve and inspect transcripts in
        **Transcript & Isoform Analysis** first.

        Then transfer the selected sequences to this module.

        Alternatively, retrieve sequences directly from NCBI
        or upload your own FASTA file.
        """
    )

    # --------------------------------------------------------
    # INPUT SOURCES
    # --------------------------------------------------------

    tab_registry, tab_ncbi, tab_fasta = st.tabs([

        "Shared Transcript Registry",

        "NCBI Retrieval",

        "FASTA Upload",

    ])

    with tab_registry:

        render_shared_registry()

    with tab_ncbi:

        render_ncbi_search()

    with tab_fasta:

        render_fasta_upload()

    # --------------------------------------------------------
    # ACTIVE RECORDS
    # --------------------------------------------------------

    records = normalize_records(

        st.session_state.primer_transcripts

    )

    st.session_state.primer_transcripts = records

    st.divider()

    if not records:

        st.warning(
            "No active transcript sequences selected."
        )

        return

    st.caption(

        "Input source: "
        + str(
            st.session_state.primer_source_mode
        )

    )

    render_active_transcripts(records)

    # --------------------------------------------------------
    # CLEAR TARGETS
    # --------------------------------------------------------

    if st.button(

        "Clear active primer targets",

        key="primer_clear_active",

    ):

        st.session_state.primer_transcripts = []

        st.session_state.primer_source_mode = ""

        st.rerun()

    # --------------------------------------------------------
    # PRIMER DESIGN
    # --------------------------------------------------------

    st.divider()

    st.subheader(
        "Primer Design and Evaluation"
    )

    st.markdown(
        """
        The RT-qPCR module performs primer generation
        and evaluation using the active transcript sequences.

        The available analysis may include:

        - Primer3-based primer design
        - Primer melting temperature
        - GC content
        - Amplicon length
        - Hairpin evaluation
        - Self-dimer evaluation
        - Primer-pair evaluation
        - Screening against the supplied transcript set
        """
    )

    st.warning(
        """
        Important:

        Screening against the supplied transcript set
        is not equivalent to genome-wide specificity analysis.

        Experimental validation remains necessary before
        using designed primers in laboratory experiments.
        """
    )

    try:

        render_rtqpcr_primer_design(records)

    except Exception as error:

        st.error(
            f"Primer design module error: {error}"
        )

        st.exception(error)


# ============================================================
# BACKWARD COMPATIBILITY
# ============================================================

def render_page():

    render_primer_design_page()