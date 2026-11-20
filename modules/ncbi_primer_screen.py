"""
NCBI BLAST-based primer specificity screening.

Important:
- This module uses the public NCBI BLAST URL API.
- It is NOT an implementation of Primer-BLAST.
- Primers are searched independently, then paired by subject,
  orientation and predicted amplicon size.
- Results are computational screening evidence, not experimental validation.
"""

from __future__ import annotations

import io
import re
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests
import streamlit as st


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

DEFAULT_BLAST_SETTINGS = {
    "program": "blastn",
    "database": "core_nt",
    "word_size": 7,
    "expect": 1000,
    "hitlist_size": 500,
    "max_poll_minutes": 15,
    "poll_interval_seconds": 60,
    "min_query_coverage": 100.0,
    "max_mismatches": 2,
    "min_amplicon": 50,
    "max_amplicon": 1000,
}


class NCBIPrimerBlastError(Exception):
    """Raised when NCBI BLAST screening fails."""


def normalize_primer(sequence: str) -> str:
    """Normalize and validate a primer sequence."""

    sequence = re.sub(r"\s+", "", str(sequence)).upper()

    if not sequence:
        raise ValueError("Primer sequence is empty.")

    if not re.fullmatch(r"[ACGT]+", sequence):
        raise ValueError(
            "Primer sequence must contain only A, C, G and T."
        )

    if len(sequence) < 15:
        raise ValueError(
            "Primer is too short for reliable BLAST screening."
        )

    return sequence


def reverse_complement(sequence: str) -> str:
    """Return reverse complement of a DNA sequence."""

    table = str.maketrans(
        "ACGT",
        "TGCA"
    )

    return sequence.translate(table)[::-1]


def build_primer_fasta(
    forward_primer: str,
    reverse_primer: str
) -> str:
    """Create a two-record FASTA query."""

    return (
        f">FORWARD_PRIMER\n"
        f"{forward_primer}\n"
        f">REVERSE_PRIMER\n"
        f"{reverse_primer}\n"
    )


def submit_blast_search(
    forward_primer: str,
    reverse_primer: str,
    email: str,
    organism: str,
    api_key: str = "",
    database: str = "core_nt",
) -> str:
    """
    Submit both primers in one BLAST request.

    Returns the NCBI Request ID (RID).
    """

    forward_primer = normalize_primer(forward_primer)
    reverse_primer = normalize_primer(reverse_primer)

    if not email or "@" not in email:
        raise ValueError(
            "A valid email address is required for NCBI requests."
        )

    query = build_primer_fasta(
        forward_primer,
        reverse_primer
    )

    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": database,
        "QUERY": query,
        "WORD_SIZE": 7,
        "EXPECT": 1000,
        "FILTER": "F",
        "HITLIST_SIZE": 500,
        "SHORT_QUERY_ADJUST": "true",
        "FORMAT_TYPE": "Tabular",
        "FORMAT_OBJECT": "Alignment",
        "ALIGNMENT_VIEW": "Tabular",
        "EMAIL": email,
        "TOOL": "BioinformaticsDashboard",
        "ENTREZ_QUERY": f"{organism}[Organism]",
    }

    if api_key:
        params["API_KEY"] = api_key

    try:
        response = requests.post(
            BLAST_URL,
            data=params,
            timeout=90
        )

        response.raise_for_status()

    except requests.RequestException as exc:
        raise NCBIPrimerBlastError(
            f"Could not submit BLAST request: {exc}"
        )

    rid_match = re.search(
        r"RID\s*=\s*([A-Z0-9-]+)",
        response.text
    )

    if not rid_match:
        raise NCBIPrimerBlastError(
            "NCBI did not return a valid RID.\n"
            + response.text[:1000]
        )

    return rid_match.group(1)


def get_blast_status(
    rid: str,
    email: str,
    api_key: str = ""
) -> Tuple[str, str]:

    params = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "Text",
        "EMAIL": email,
        "TOOL": "BioinformaticsDashboard",
    }

    if api_key:
        params["API_KEY"] = api_key

    try:
        response = requests.get(
            BLAST_URL,
            params=params,
            timeout=90
        )

        response.raise_for_status()

    except requests.RequestException as exc:
        raise NCBIPrimerBlastError(
            f"Could not retrieve BLAST status: {exc}"
        )

    text = response.text

    if "Status=WAITING" in text:
        return "WAITING", text

    if "Status=FAILED" in text:
        return "FAILED", text

    if "Status=UNKNOWN" in text:
        return "UNKNOWN", text

    if "Status=READY" in text:
        return "READY", text

    # Some responses may directly contain result text.
    if "Sequences producing significant alignments" in text:
        return "READY", text

    return "UNKNOWN", text


def wait_for_blast_results(
    rid: str,
    email: str,
    api_key: str = "",
    max_poll_minutes: int = 15,
    poll_interval_seconds: int = 60,
    progress_callback=None,
) -> str:
    """
    Wait for BLAST completion.

    NCBI recommends avoiding frequent polling.
    This implementation uses a 60-second default interval.
    """

    start_time = time.time()
    max_wait = max_poll_minutes * 60

    while time.time() - start_time < max_wait:

        status, response_text = get_blast_status(
            rid=rid,
            email=email,
            api_key=api_key
        )

        elapsed = time.time() - start_time

        if progress_callback:
            progress_callback(
                min(elapsed / max_wait, 1.0),
                status
            )

        if status == "READY":
            return response_text

        if status in ("FAILED", "UNKNOWN"):
            raise NCBIPrimerBlastError(
                f"BLAST job status: {status}\n"
                + response_text[:1000]
            )

        time.sleep(poll_interval_seconds)

    raise NCBIPrimerBlastError(
        "BLAST job timed out. "
        "You can retry the screening later."
    )


def parse_blast_tabular(
    blast_text: str
) -> pd.DataFrame:
    """
    Parse standard 12-column BLAST tabular output.

    Expected columns:
    qseqid, sseqid, pident, length, mismatch, gapopen,
    qstart, qend, sstart, send, evalue, bitscore
    """

    columns = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]

    lines = []

    for line in blast_text.splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        fields = line.split("\t")

        if len(fields) >= 12:
            lines.append(fields[:12])

    if not lines:
        return pd.DataFrame(columns=columns)

    df = pd.DataFrame(
        lines,
        columns=columns
    )

    numeric_columns = [
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce"
        )

    df = df.dropna(
        subset=[
            "pident",
            "length",
            "sstart",
            "send",
        ]
    )

    return df


def filter_full_length_hits(
    hits: pd.DataFrame,
    primer_length: int,
    max_mismatches: int = 2,
    min_query_coverage: float = 100.0,
) -> pd.DataFrame:
    """
    Retain sufficiently complete primer alignments.

    Note:
    BLAST local alignment may omit terminal mismatches.
    Therefore, query coverage and mismatch count are both checked.
    """

    if hits.empty:
        return hits.copy()

    result = hits.copy()

    result["query_coverage"] = (
        result["length"] / primer_length
    ) * 100

    result = result[
        (result["query_coverage"] >= min_query_coverage)
        &
        (result["mismatch"] <= max_mismatches)
        &
        (result["gapopen"] == 0)
    ]

    return result.copy()


def extract_primer_hits(
    hits: pd.DataFrame,
    query_name: str
) -> pd.DataFrame:

    if hits.empty:
        return hits.copy()

    result = hits[
        hits["qseqid"].str.contains(
            query_name,
            case=False,
            na=False
        )
    ].copy()

    return result


def predict_amplicons(
    forward_hits: pd.DataFrame,
    reverse_hits: pd.DataFrame,
    min_amplicon: int = 50,
    max_amplicon: int = 1000,
) -> pd.DataFrame:
    """
    Pair forward and reverse primer hits.

    A candidate amplicon requires:
    - Same subject sequence
    - Forward primer on the plus strand
    - Reverse primer on the minus strand
    - Forward binding site upstream of reverse binding site
    - Product length within requested range
    """

    results = []

    if forward_hits.empty or reverse_hits.empty:
        return pd.DataFrame()

    for _, f_hit in forward_hits.iterrows():

        for _, r_hit in reverse_hits.iterrows():

            if f_hit["sseqid"] != r_hit["sseqid"]:
                continue

            # Forward primer must align in plus orientation.
            forward_plus = (
                f_hit["sstart"] < f_hit["send"]
            )

            # Reverse primer must align in minus orientation.
            reverse_minus = (
                r_hit["sstart"] > r_hit["send"]
            )

            if not forward_plus or not reverse_minus:
                continue

            f_start = int(
                min(f_hit["sstart"], f_hit["send"])
            )

            f_end = int(
                max(f_hit["sstart"], f_hit["send"])
            )

            r_start = int(
                min(r_hit["sstart"], r_hit["send"])
            )

            r_end = int(
                max(r_hit["sstart"], r_hit["send"])
            )

            if r_start <= f_start:
                continue

            product_size = r_end - f_start + 1

            if not (
                min_amplicon
                <= product_size
                <= max_amplicon
            ):
                continue

            results.append({
                "subject_accession": f_hit["sseqid"],
                "forward_start": f_start,
                "forward_end": f_end,
                "reverse_start": r_start,
                "reverse_end": r_end,
                "amplicon_size": product_size,
                "forward_mismatches": int(
                    f_hit["mismatch"]
                ),
                "reverse_mismatches": int(
                    r_hit["mismatch"]
                ),
                "forward_identity": float(
                    f_hit["pident"]
                ),
                "reverse_identity": float(
                    r_hit["pident"]
                ),
                "forward_bitscore": float(
                    f_hit["bitscore"]
                ),
                "reverse_bitscore": float(
                    r_hit["bitscore"]
                ),
            })

    return pd.DataFrame(results)


def run_primer_specificity_screen(
    forward_primer: str,
    reverse_primer: str,
    email: str,
    organism: str,
    api_key: str = "",
    database: str = "core_nt",
    progress_callback=None,
) -> Dict:

    forward_primer = normalize_primer(
        forward_primer
    )

    reverse_primer = normalize_primer(
        reverse_primer
    )

    rid = submit_blast_search(
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
        email=email,
        organism=organism,
        api_key=api_key,
        database=database,
    )

    blast_text = wait_for_blast_results(
        rid=rid,
        email=email,
        api_key=api_key,
        progress_callback=progress_callback,
    )

    all_hits = parse_blast_tabular(
        blast_text
    )

    forward_hits = extract_primer_hits(
        all_hits,
        "FORWARD_PRIMER"
    )

    reverse_hits = extract_primer_hits(
        all_hits,
        "REVERSE_PRIMER"
    )

    forward_hits = filter_full_length_hits(
        forward_hits,
        primer_length=len(forward_primer),
        max_mismatches=2,
        min_query_coverage=100,
    )

    reverse_hits = filter_full_length_hits(
        reverse_hits,
        primer_length=len(reverse_primer),
        max_mismatches=2,
        min_query_coverage=100,
    )

    amplicons = predict_amplicons(
        forward_hits,
        reverse_hits,
        min_amplicon=50,
        max_amplicon=1000,
    )

    return {
        "rid": rid,
        "all_hits": all_hits,
        "forward_hits": forward_hits,
        "reverse_hits": reverse_hits,
        "amplicons": amplicons,
        "blast_raw": blast_text,
    }


def render_ncbi_specificity_screen(
    forward_primer: str,
    reverse_primer: str,
    organism: str,
):
    """
    Streamlit UI for NCBI primer specificity screening.
    """

    st.subheader(
        "NCBI BLAST — Primer Specificity Screening"
    )

    st.warning(
        "This section uses the public NCBI BLAST and is not equivalent to Primer-BLAST. Its results represent a computational screen, not a definitive confirmation of specificity. "       
    )

    email = st.text_input(
        "NCBI contact email",
        key="ncbi_primer_email",
        help="An email address is required for responsible NCBI API use."
    )

    api_key = st.text_input(
        "NCBI API key (optional)",
        type="password",
        key="ncbi_primer_api_key"
    )

    database = st.selectbox(
        "BLAST database",
        options=[
            "core_nt",
            "refseq_rna",
        ],
        index=0,
        key="ncbi_primer_database"
    )

    if st.button(
        "Run NCBI Specificity Screening",
        type="primary",
        key="run_ncbi_primer_screen"
    ):

        if not email or "@" not in email:
            st.error(
                "Please enter a valid email address."
            )
            return

        try:

            progress = st.progress(0)
            status_box = st.empty()

            def update_progress(value, status):

                progress.progress(
                    min(max(value, 0.0), 1.0)
                )

                status_box.info(
                    f"NCBI BLAST status: {status}"
                )

            with st.spinner(
                "Submitting primers and waiting for NCBI..."
            ):

                result = run_primer_specificity_screen(
                    forward_primer=forward_primer,
                    reverse_primer=reverse_primer,
                    email=email,
                    organism=organism,
                    api_key=api_key,
                    database=database,
                    progress_callback=update_progress,
                )

            st.session_state[
                "ncbi_primer_screen_result"
            ] = result

            progress.progress(1.0)

            st.success(
                f"BLAST completed. RID: {result['rid']}"
            )

        except Exception as exc:

            st.error(
                f"Specificity screening failed: {exc}"
            )

    result = st.session_state.get(
        "ncbi_primer_screen_result"
    )

    if not result:
        return

    st.divider()

    st.markdown("### Screening summary")

    col1, col2, col3 = st.columns(3)

    col1.metric(
        "Forward hits",
        len(result["forward_hits"])
    )

    col2.metric(
        "Reverse hits",
        len(result["reverse_hits"])
    )

    col3.metric(
        "Predicted amplicons",
        len(result["amplicons"])
    )

    if result["amplicons"].empty:

        st.info(
            "No paired amplicons were detected under "
            "the current alignment and product-size criteria. "
            "This does not prove that the primers are specific."
        )

    else:

        st.markdown(
            "### Predicted paired amplicons"
        )

        st.dataframe(
            result["amplicons"],
            use_container_width=True,
            hide_index=True
        )

        csv = result["amplicons"].to_csv(
            index=False
        ).encode("utf-8")

        st.download_button(
            "Download amplicon results",
            data=csv,
            file_name="ncbi_primer_amplicons.csv",
            mime="text/csv"
        )

    with st.expander(
        "Forward primer BLAST hits"
    ):

        st.dataframe(
            result["forward_hits"],
            use_container_width=True
        )

    with st.expander(
        "Reverse primer BLAST hits"
    ):

        st.dataframe(
            result["reverse_hits"],
            use_container_width=True
        )

    with st.expander(
        "Raw BLAST output"
    ):

        st.text(
            result["blast_raw"][:20000]
        )