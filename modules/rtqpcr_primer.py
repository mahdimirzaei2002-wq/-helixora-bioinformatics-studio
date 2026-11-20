import re
import primer3
import pandas as pd
import streamlit as st

from Bio.Seq import Seq

from modules.ncbi_primer_screen import (
    render_ncbi_specificity_screen
)


# =========================================================
# 1. DEFAULT PRIMER DESIGN PARAMETERS
# =========================================================

DEFAULT_SETTINGS = {

    "PRIMER_NUM_RETURN": 10,

    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MAX_SIZE": 24,

    "PRIMER_MIN_TM": 58.0,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MAX_TM": 62.0,

    "PRIMER_PAIR_MAX_DIFF_TM": 2.0,

    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,

    "PRIMER_PRODUCT_SIZE_RANGE": [[70, 200]],

    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_GC_CLAMP": 1,

    "PRIMER_MAX_SELF_ANY": 8.0,
    "PRIMER_MAX_SELF_END": 3.0,

    "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
    "PRIMER_PAIR_MAX_COMPL_END": 3.0,

    "PRIMER_EXPLAIN_FLAG": 1,
}


# =========================================================
# 2. TRANSCRIPT NORMALIZATION
# =========================================================

def normalize_transcript(item, index=0):

    accession = f"Transcript_{index + 1}"
    description = ""

    if hasattr(item, "seq"):

        sequence = str(item.seq)

        accession = (
            getattr(item, "id", None)
            or accession
        )

        description = getattr(
            item,
            "description",
            ""
        )

    elif isinstance(item, dict):

        accession = (
            item.get("accession")
            or item.get("id")
            or item.get("label")
            or accession
        )

        description = item.get(
            "title",
            ""
        )

        record = item.get("record")

        if hasattr(record, "seq"):

            sequence = str(record.seq)

        else:

            sequence = (
                item.get("sequence")
                or item.get("seq")
                or ""
            )

    else:

        raise ValueError(
            f"Unsupported transcript format: {type(item)}"
        )

    sequence = re.sub(
        r"\s+",
        "",
        sequence
    ).upper()

    sequence = sequence.replace("U", "T")

    if not sequence:

        raise ValueError(
            f"No sequence found for {accession}"
        )

    return {

        "accession": str(accession),

        "description": str(description),

        "sequence": sequence

    }


# =========================================================
# 3. TEMPLATE VALIDATION
# =========================================================

def validate_template(sequence):

    sequence = sequence.upper()

    if len(sequence) < 100:

        raise ValueError(
            "Transcript is too short for RT-qPCR primer design."
        )

    invalid = set(sequence) - set("ACGT")

    if invalid:

        raise ValueError(
            "The selected transcript contains ambiguous "
            f"nucleotides: {', '.join(sorted(invalid))}. "
            "Select a transcript containing only A, C, G and T."
        )

    return True


# =========================================================
# 4. PRIMER DESIGN USING PRIMER3
# =========================================================

def design_primer_pairs(
    sequence,
    sequence_id="Target_transcript",
    settings=None
):

    validate_template(sequence)

    parameters = DEFAULT_SETTINGS.copy()

    if settings:

        parameters.update(settings)

    sequence_args = {

        "SEQUENCE_ID": sequence_id,

        "SEQUENCE_TEMPLATE": sequence,

    }

    result = primer3.bindings.design_primers(
        sequence_args,
        parameters
    )

    return result


# =========================================================
# 5. GC CONTENT
# =========================================================

def calculate_gc(sequence):

    sequence = sequence.upper()

    if not sequence:

        return 0.0

    gc_count = (
        sequence.count("G")
        + sequence.count("C")
    )

    return (
        gc_count / len(sequence)
    ) * 100


# =========================================================
# 6. SAFE SECONDARY STRUCTURE ANALYSIS
# =========================================================

def safe_structure_analysis(
    sequence,
    structure_function
):

    try:

        result = structure_function(sequence)

        return {

            "Tm": float(result.tm),

            "Delta G": float(result.dg),

            "Structure Found": bool(
                result.structure_found
            )

        }

    except Exception:

        return {

            "Tm": None,

            "Delta G": None,

            "Structure Found": None

        }


# =========================================================
# 7. PRIMER PAIR THERMODYNAMIC EVALUATION
# =========================================================

def evaluate_primer_pair(
    forward_primer,
    reverse_primer,
    product_size=None
):

    forward_primer = forward_primer.upper()

    reverse_primer = reverse_primer.upper()

    # -----------------------------------------------------
    # Basic properties
    # -----------------------------------------------------

    f_length = len(forward_primer)

    r_length = len(reverse_primer)

    f_gc = calculate_gc(forward_primer)

    r_gc = calculate_gc(reverse_primer)

    # -----------------------------------------------------
    # Melting temperature
    # -----------------------------------------------------

    try:

        f_tm = float(
            primer3.calc_tm(forward_primer)
        )

        r_tm = float(
            primer3.calc_tm(reverse_primer)
        )

        tm_difference = abs(
            f_tm - r_tm
        )

    except Exception:

        f_tm = None

        r_tm = None

        tm_difference = None

    # -----------------------------------------------------
    # Hairpin
    # -----------------------------------------------------

    f_hairpin = safe_structure_analysis(
        forward_primer,
        primer3.calc_hairpin
    )

    r_hairpin = safe_structure_analysis(
        reverse_primer,
        primer3.calc_hairpin
    )

    # -----------------------------------------------------
    # Homodimer
    # -----------------------------------------------------

    f_homodimer = safe_structure_analysis(
        forward_primer,
        primer3.calc_homodimer
    )

    r_homodimer = safe_structure_analysis(
        reverse_primer,
        primer3.calc_homodimer
    )

    # -----------------------------------------------------
    # Heterodimer
    # -----------------------------------------------------

    heterodimer = {

        "Tm": None,

        "Delta G": None,

        "Structure Found": None

    }

    try:

        heterodimer_result = primer3.calc_heterodimer(
            forward_primer,
            reverse_primer
        )

        heterodimer = {

            "Tm": float(
                heterodimer_result.tm
            ),

            "Delta G": float(
                heterodimer_result.dg
            ),

            "Structure Found": bool(
                heterodimer_result.structure_found
            )

        }

    except Exception:

        pass

    # -----------------------------------------------------
    # Build metrics
    # -----------------------------------------------------

    metrics = {

        "Forward Length": f_length,

        "Reverse Length": r_length,

        "Forward GC (%)": round(f_gc, 2),

        "Reverse GC (%)": round(r_gc, 2),

        "Forward Tm (°C)": (
            round(f_tm, 2)
            if f_tm is not None
            else None
        ),

        "Reverse Tm (°C)": (
            round(r_tm, 2)
            if r_tm is not None
            else None
        ),

        "Tm Difference (°C)": (
            round(tm_difference, 2)
            if tm_difference is not None
            else None
        ),

        "Forward Hairpin Tm": f_hairpin["Tm"],

        "Reverse Hairpin Tm": r_hairpin["Tm"],

        "Forward Homodimer Tm": f_homodimer["Tm"],

        "Reverse Homodimer Tm": r_homodimer["Tm"],

        "Heterodimer Tm": heterodimer["Tm"],

        "Forward Hairpin ΔG": f_hairpin["Delta G"],

        "Reverse Hairpin ΔG": r_hairpin["Delta G"],

        "Forward Homodimer ΔG": f_homodimer["Delta G"],

        "Reverse Homodimer ΔG": r_homodimer["Delta G"],

        "Heterodimer ΔG": heterodimer["Delta G"],

        "Product Size": product_size,

    }

    return metrics


# =========================================================
# 8. PRIMER QUALITY SCORE
# =========================================================

def calculate_primer_quality_score(metrics):

    score = 100.0

    warnings = []

    # -----------------------------------------------------
    # Tm difference
    # -----------------------------------------------------

    tm_diff = metrics.get(
        "Tm Difference (°C)"
    )

    if tm_diff is None:

        score -= 20

        warnings.append(
            "Tm could not be calculated."
        )

    elif tm_diff > 3:

        score -= 25

        warnings.append(
            "Large difference between primer Tm values."
        )

    elif tm_diff > 2:

        score -= 15

        warnings.append(
            "Primer Tm difference is above 2°C."
        )

    elif tm_diff > 1:

        score -= 5

    # -----------------------------------------------------
    # GC content
    # -----------------------------------------------------

    for key in [
        "Forward GC (%)",
        "Reverse GC (%)"
    ]:

        gc = metrics.get(key)

        if gc is None:

            score -= 10

            warnings.append(
                f"{key} could not be evaluated."
            )

        elif gc < 40 or gc > 60:

            score -= 15

            warnings.append(
                f"{key} is outside the preferred 40–60% range."
            )

        elif gc < 45 or gc > 55:

            score -= 4

    # -----------------------------------------------------
    # Primer length
    # -----------------------------------------------------

    for key in [
        "Forward Length",
        "Reverse Length"
    ]:

        length = metrics.get(key)

        if length is None:

            score -= 10

        elif length < 18 or length > 24:

            score -= 10

            warnings.append(
                f"{key} is outside the preferred 18–24 nt range."
            )

    # -----------------------------------------------------
    # Hairpin
    # -----------------------------------------------------

    for key in [
        "Forward Hairpin Tm",
        "Reverse Hairpin Tm"
    ]:

        tm = metrics.get(key)

        if tm is not None:

            if tm >= 50:

                score -= 15

                warnings.append(
                    f"Potentially strong hairpin: {key}."
                )

            elif tm >= 40:

                score -= 7

                warnings.append(
                    f"Potential hairpin concern: {key}."
                )

    # -----------------------------------------------------
    # Homodimer
    # -----------------------------------------------------

    for key in [
        "Forward Homodimer Tm",
        "Reverse Homodimer Tm"
    ]:

        tm = metrics.get(key)

        if tm is not None:

            if tm >= 50:

                score -= 15

                warnings.append(
                    f"Potentially strong homodimer: {key}."
                )

            elif tm >= 40:

                score -= 7

                warnings.append(
                    f"Potential homodimer concern: {key}."
                )

    # -----------------------------------------------------
    # Heterodimer
    # -----------------------------------------------------

    heterodimer_tm = metrics.get(
        "Heterodimer Tm"
    )

    if heterodimer_tm is not None:

        if heterodimer_tm >= 50:

            score -= 20

            warnings.append(
                "Potentially strong forward-reverse heterodimer."
            )

        elif heterodimer_tm >= 40:

            score -= 10

            warnings.append(
                "Potential heterodimer concern."
            )

    score = max(
        0,
        min(100, score)
    )

    return round(score, 2), warnings


# =========================================================
# 9. TRANSCRIPT ISOFORM SCREENING
# =========================================================

def screen_transcripts(
    forward,
    reverse,
    transcripts,
    min_product=70,
    max_product=200
):

    """
    Exact-match screening against retrieved transcripts.

    This is not genome-wide specificity analysis.
    """

    forward = forward.upper()

    reverse = reverse.upper()

    reverse_binding = str(
        Seq(reverse).reverse_complement()
    )

    results = []

    for transcript in transcripts:

        sequence = transcript["sequence"]

        forward_positions = [

            match.start()

            for match in re.finditer(
                f"(?={re.escape(forward)})",
                sequence
            )

        ]

        reverse_positions = [

            match.start()

            for match in re.finditer(
                f"(?={re.escape(reverse_binding)})",
                sequence
            )

        ]

        products = []

        for f_pos in forward_positions:

            for r_pos in reverse_positions:

                if r_pos <= f_pos:

                    continue

                product_size = (
                    r_pos
                    + len(reverse_binding)
                    - f_pos
                )

                if min_product <= product_size <= max_product:

                    products.append(product_size)

        results.append({

            "Accession": transcript["accession"],

            "Exact Forward Matches": len(
                forward_positions
            ),

            "Exact Reverse Matches": len(
                reverse_positions
            ),

            "Amplicon Found": (
                "Yes"
                if products
                else "No"
            ),

            "Amplicon Sizes (bp)": (
                ", ".join(
                    map(
                        str,
                        sorted(set(products))
                    )
                )
                if products
                else "-"
            ),

            "Number of Amplicons": len(products),

        })

    return pd.DataFrame(results)


# =========================================================
# 10. EXTRACT PRIMER3 CANDIDATES
# =========================================================

def extract_primer_candidates(
    result,
    transcripts,
    min_product=70,
    max_product=200
):

    candidates = []

    number = result.get(
        "PRIMER_PAIR_NUM_RETURNED",
        0
    )

    for i in range(number):

        forward = result.get(
            f"PRIMER_LEFT_{i}_SEQUENCE"
        )

        reverse = result.get(
            f"PRIMER_RIGHT_{i}_SEQUENCE"
        )

        if not forward or not reverse:

            continue

        forward_tm = result.get(
            f"PRIMER_LEFT_{i}_TM"
        )

        reverse_tm = result.get(
            f"PRIMER_RIGHT_{i}_TM"
        )

        forward_gc = result.get(
            f"PRIMER_LEFT_{i}_GC_PERCENT"
        )

        reverse_gc = result.get(
            f"PRIMER_RIGHT_{i}_GC_PERCENT"
        )

        product_size = result.get(
            f"PRIMER_PAIR_{i}_PRODUCT_SIZE"
        )

        penalty = result.get(
            f"PRIMER_PAIR_{i}_PENALTY"
        )

        screening = screen_transcripts(
            forward,
            reverse,
            transcripts,
            min_product,
            max_product
        )

        matched = screening.loc[
            screening["Amplicon Found"] == "Yes",
            "Accession"
        ].tolist()

        candidates.append({

            "Candidate": i + 1,

            "Forward Primer": forward,

            "Reverse Primer": reverse,

            "Forward Length": len(forward),

            "Reverse Length": len(reverse),

            "Forward Tm (°C)": round(
                float(forward_tm), 2
            ),

            "Reverse Tm (°C)": round(
                float(reverse_tm), 2
            ),

            "Tm Difference (°C)": round(
                abs(
                    float(forward_tm)
                    - float(reverse_tm)
                ),
                2
            ),

            "Forward GC (%)": round(
                float(forward_gc), 2
            ),

            "Reverse GC (%)": round(
                float(reverse_gc), 2
            ),

            "Amplicon Size (bp)": int(
                product_size
            ),

            "Primer3 Penalty": round(
                float(penalty), 3
            ),

            "Transcript Matches": ", ".join(matched),

            "Matched Transcript Count": len(matched)

        })

    return candidates


# =========================================================
# 11. ENRICH AND RANK CANDIDATES
# =========================================================

def enrich_and_rank_candidates(candidates):

    evaluated_candidates = []

    for candidate in candidates:

        candidate = candidate.copy()

        forward = candidate.get(
            "Forward Primer"
        )

        reverse = candidate.get(
            "Reverse Primer"
        )

        if not forward or not reverse:

            continue

        product_size = candidate.get(
            "Amplicon Size (bp)"
        )

        metrics = evaluate_primer_pair(
            forward,
            reverse,
            product_size
        )

        score, warnings = calculate_primer_quality_score(
            metrics
        )

        candidate.update(metrics)

        candidate["Quality Score"] = score

        candidate["Warnings"] = warnings

        evaluated_candidates.append(candidate)

    evaluated_candidates.sort(
        key=lambda x: (
            x["Quality Score"],
            -x.get("Primer3 Penalty", 0)
        ),
        reverse=True
    )

    for index, candidate in enumerate(
        evaluated_candidates,
        start=1
    ):

        candidate["Design Rank"] = index

    return evaluated_candidates


# =========================================================
# 12. STREAMLIT INTERFACE
# =========================================================

def render_rtqpcr_primer_design(transcript_items):

    st.header("RT-qPCR Primer Design")

    st.write(
        "Design, evaluate and screen candidate RT-qPCR "
        "primer pairs using Primer3 and NCBI BLAST."
    )

    # =====================================================
    # STEP 1 — VALIDATE TRANSCRIPTS
    # =====================================================

    if not transcript_items:

        st.warning(
            "No retrieved transcript records were found. "
            "First search for a gene and retrieve its transcripts."
        )

        return

    normalized = []

    for index, item in enumerate(transcript_items):

        try:

            normalized.append(
                normalize_transcript(item, index)
            )

        except Exception as error:

            st.warning(
                f"Transcript {index + 1} skipped: {error}"
            )

    if not normalized:

        st.error(
            "No valid transcript sequences are available."
        )

        return

    # =====================================================
    # STEP 2 — SELECT REFERENCE TRANSCRIPT
    # =====================================================

    st.subheader("1. Reference Transcript Selection")

    transcript_labels = [

        f"{t['accession']} | {len(t['sequence'])} bp"

        for t in normalized

    ]

    selected_index = st.selectbox(
        "Select reference transcript",
        range(len(normalized)),
        format_func=lambda x: transcript_labels[x],
        key="rtqpcr_reference_transcript"
    )

    selected = normalized[selected_index]

    st.info(
        f"Selected transcript: {selected['accession']} "
        f"({len(selected['sequence'])} bp)"
    )

    # =====================================================
    # STEP 3 — DESIGN PARAMETERS
    # =====================================================

    st.divider()

    st.subheader("2. Primer Design Parameters")

    col1, col2, col3 = st.columns(3)

    with col1:

        min_size = st.number_input(
            "Minimum primer length",
            min_value=15,
            max_value=30,
            value=18,
            key="rtqpcr_min_size"
        )

        opt_size = st.number_input(
            "Optimal primer length",
            min_value=15,
            max_value=30,
            value=20,
            key="rtqpcr_opt_size"
        )

        max_size = st.number_input(
            "Maximum primer length",
            min_value=15,
            max_value=35,
            value=24,
            key="rtqpcr_max_size"
        )

    with col2:

        min_tm = st.number_input(
            "Minimum Tm (°C)",
            value=58.0,
            step=0.5,
            key="rtqpcr_min_tm"
        )

        opt_tm = st.number_input(
            "Optimal Tm (°C)",
            value=60.0,
            step=0.5,
            key="rtqpcr_opt_tm"
        )

        max_tm = st.number_input(
            "Maximum Tm (°C)",
            value=62.0,
            step=0.5,
            key="rtqpcr_max_tm"
        )

    with col3:

        min_gc = st.number_input(
            "Minimum GC (%)",
            value=40.0,
            step=1.0,
            key="rtqpcr_min_gc"
        )

        max_gc = st.number_input(
            "Maximum GC (%)",
            value=60.0,
            step=1.0,
            key="rtqpcr_max_gc"
        )

        number_return = st.number_input(
            "Number of candidates",
            min_value=1,
            max_value=20,
            value=10,
            key="rtqpcr_number_return"
        )

    col4, col5 = st.columns(2)

    with col4:

        min_product = st.number_input(
            "Minimum amplicon size (bp)",
            min_value=50,
            max_value=500,
            value=70,
            key="rtqpcr_min_product"
        )

    with col5:

        max_product = st.number_input(
            "Maximum amplicon size (bp)",
            min_value=50,
            max_value=1000,
            value=200,
            key="rtqpcr_max_product"
        )

    # =====================================================
    # STEP 4 — VALIDATE PARAMETERS
    # =====================================================

    valid_parameters = (

        min_size <= opt_size <= max_size

        and min_tm <= opt_tm <= max_tm

        and min_gc < max_gc

        and min_product < max_product

    )

    if not valid_parameters:

        st.error(
            "Invalid parameter ranges. "
            "Please check minimum, optimal and maximum values."
        )

        return

    # =====================================================
    # STEP 5 — BUILD PRIMER3 SETTINGS
    # =====================================================

    settings = {

        "PRIMER_NUM_RETURN": int(number_return),

        "PRIMER_MIN_SIZE": int(min_size),

        "PRIMER_OPT_SIZE": int(opt_size),

        "PRIMER_MAX_SIZE": int(max_size),

        "PRIMER_MIN_TM": float(min_tm),

        "PRIMER_OPT_TM": float(opt_tm),

        "PRIMER_MAX_TM": float(max_tm),

        "PRIMER_PAIR_MAX_DIFF_TM": 2.0,

        "PRIMER_MIN_GC": float(min_gc),

        "PRIMER_MAX_GC": float(max_gc),

        "PRIMER_PRODUCT_SIZE_RANGE": [
            [int(min_product), int(max_product)]
        ],

    }

    # =====================================================
    # STEP 6 — RUN PRIMER3
    # =====================================================

    st.divider()

    st.subheader("3. Primer Design")

    if st.button(
        "Design RT-qPCR Primers",
        type="primary",
        key="run_rtqpcr_design"
    ):

        try:

            with st.spinner(
                "Primer3 is designing candidate primer pairs..."
            ):

                result = design_primer_pairs(
                    selected["sequence"],
                    selected["accession"],
                    settings
                )

                candidates = extract_primer_candidates(
                    result,
                    normalized,
                    min_product,
                    max_product
                )

                candidates = enrich_and_rank_candidates(
                    candidates
                )

            if not candidates:

                st.session_state["rtqpcr_candidates"] = []

                st.session_state["rtqpcr_screening"] = {}

                st.session_state["rtqpcr_selected_accession"] = (
                    selected["accession"]
                )

                st.error(
                    "Primer3 could not identify suitable "
                    "primer pairs with the selected parameters."
                )

                st.write(
                    "Forward primer:",
                    result.get(
                        "PRIMER_LEFT_EXPLAIN",
                        "Not available"
                    )
                )

                st.write(
                    "Reverse primer:",
                    result.get(
                        "PRIMER_RIGHT_EXPLAIN",
                        "Not available"
                    )
                )

                st.write(
                    "Primer pairs:",
                    result.get(
                        "PRIMER_PAIR_EXPLAIN",
                        "Not available"
                    )
                )

                return

            # Save candidates

            st.session_state["rtqpcr_candidates"] = candidates

            st.session_state["rtqpcr_selected_accession"] = (
                selected["accession"]
            )

            st.session_state["rtqpcr_selected_candidate"] = (
                str(candidates[0]["Candidate"])
            )

            # Save transcript screening

            screening_results = {}

            for candidate in candidates:

                candidate_id = str(
                    candidate["Candidate"]
                )

                screening_results[candidate_id] = (
                    screen_transcripts(
                        candidate["Forward Primer"],
                        candidate["Reverse Primer"],
                        normalized,
                        min_product,
                        max_product
                    )
                )

            st.session_state["rtqpcr_screening"] = (
                screening_results
            )

            st.success(
                f"Successfully generated {len(candidates)} "
                "candidate primer pairs."
            )

        except Exception as error:

            st.error(
                f"Primer design failed: {error}"
            )

            return

    # =====================================================
    # STEP 7 — LOAD CANDIDATES
    # =====================================================

    candidates = st.session_state.get(
        "rtqpcr_candidates",
        []
    )

    saved_accession = st.session_state.get(
        "rtqpcr_selected_accession"
    )

    if saved_accession != selected["accession"]:

        st.info(
            "The reference transcript has changed. "
            "Please run primer design again."
        )

        return

    if not candidates:

        st.info(
            "Run Primer3 to generate candidate primer pairs."
        )

        return

    # =====================================================
    # STEP 8 — DISPLAY CANDIDATES
    # =====================================================

    st.divider()

    st.subheader("4. Candidate Primer Pairs")

    results_df = pd.DataFrame(candidates)

    display_columns = [

        "Design Rank",

        "Candidate",

        "Forward Primer",

        "Reverse Primer",

        "Forward Tm (°C)",

        "Reverse Tm (°C)",

        "Tm Difference (°C)",

        "Forward GC (%)",

        "Reverse GC (%)",

        "Amplicon Size (bp)",

        "Quality Score",

        "Matched Transcript Count"

    ]

    available_columns = [

        column

        for column in display_columns

        if column in results_df.columns

    ]

    st.dataframe(
        results_df[available_columns],
        use_container_width=True,
        hide_index=True
    )

    # =====================================================
    # STEP 9 — SELECT CANDIDATE
    # =====================================================

    st.divider()

    st.subheader("5. Detailed Primer Evaluation")

    candidate_ids = [

        str(candidate["Candidate"])

        for candidate in candidates

    ]

    current_candidate = st.session_state.get(
        "rtqpcr_selected_candidate",
        candidate_ids[0]
    )

    if current_candidate not in candidate_ids:

        current_candidate = candidate_ids[0]

    selected_candidate_id = st.selectbox(
        "Select a primer pair",
        candidate_ids,
        index=candidate_ids.index(current_candidate),
        key="rtqpcr_selected_candidate"
    )

    selected_candidate = next(
        candidate
        for candidate in candidates
        if str(candidate["Candidate"]) == selected_candidate_id
    )

    forward_primer = selected_candidate["Forward Primer"]

    reverse_primer = selected_candidate["Reverse Primer"]

    # =====================================================
    # STEP 10 — SELECTED PRIMERS
    # =====================================================

    st.markdown("### Selected Primer Pair")

    primer_col1, primer_col2 = st.columns(2)

    with primer_col1:

        st.markdown("**Forward Primer**")

        st.code(
            forward_primer,
            language=None
        )

    with primer_col2:

        st.markdown("**Reverse Primer**")

        st.code(
            reverse_primer,
            language=None
        )

    # =====================================================
    # STEP 11 — QUALITY SCORE
    # =====================================================

    st.divider()

    st.subheader("6. Primer Pair Quality Assessment")

    quality_score = selected_candidate.get(
        "Quality Score",
        0
    )

    rank = selected_candidate.get(
        "Design Rank",
        "-"
    )

    col1, col2 = st.columns(2)

    col1.metric(
        "Quality Score",
        f"{quality_score}/100"
    )

    col2.metric(
        "Candidate Rank",
        rank
    )

    st.progress(
        int(quality_score)
    )

    st.caption(
        "The quality score is a heuristic ranking metric "
        "and does not represent the probability of experimental success."
    )

    # =====================================================
    # STEP 12 — THERMODYNAMIC PROPERTIES
    # =====================================================

    st.divider()

    st.subheader("7. Thermodynamic Properties")

    metric_keys = [

        "Forward Length",

        "Reverse Length",

        "Forward GC (%)",

        "Reverse GC (%)",

        "Forward Tm (°C)",

        "Reverse Tm (°C)",

        "Tm Difference (°C)",

    ]

    metric_data = [

        {

            "Parameter": key,

            "Value": selected_candidate.get(
                key,
                "Not available"
            )

        }

        for key in metric_keys

    ]

    st.dataframe(
        pd.DataFrame(metric_data),
        use_container_width=True,
        hide_index=True
    )

    # =====================================================
    # STEP 13 — SECONDARY STRUCTURE ANALYSIS
    # =====================================================

    st.divider()

    st.subheader("8. Secondary Structure Analysis")

    structure_keys = [

        "Forward Hairpin Tm",

        "Reverse Hairpin Tm",

        "Forward Homodimer Tm",

        "Reverse Homodimer Tm",

        "Heterodimer Tm",

        "Forward Hairpin ΔG",

        "Reverse Hairpin ΔG",

        "Forward Homodimer ΔG",

        "Reverse Homodimer ΔG",

        "Heterodimer ΔG",

    ]

    structure_data = [

        {

            "Structure Parameter": key,

            "Value": selected_candidate.get(
                key,
                "Not available"
            )

        }

        for key in structure_keys

    ]

    st.dataframe(
        pd.DataFrame(structure_data),
        use_container_width=True,
        hide_index=True
    )

    # =====================================================
    # STEP 14 — QUALITY WARNINGS
    # =====================================================

    st.divider()

    st.subheader("9. Quality Warnings")

    warnings = selected_candidate.get(
        "Warnings",
        []
    )

    if warnings:

        for warning in warnings:

            st.warning(warning)

    else:

        st.success(
            "No major concerns were detected by the "
            "current heuristic screening criteria."
        )

    # =====================================================
    # STEP 15 — TRANSCRIPT COVERAGE
    # =====================================================

    st.divider()

    st.subheader("10. Transcript Coverage Analysis")

    screening_results = st.session_state.get(
        "rtqpcr_screening",
        {}
    )

    candidate_screening = screening_results.get(
        selected_candidate_id
    )

    if candidate_screening is not None:

        st.dataframe(
            candidate_screening,
            use_container_width=True,
            hide_index=True
        )

        screening_csv = candidate_screening.to_csv(
            index=False
        ).encode("utf-8-sig")

        st.download_button(
            "Download Transcript Screening (CSV)",
            data=screening_csv,
            file_name=(
                f"transcript_screening_"
                f"candidate_{selected_candidate_id}.csv"
            ),
            mime="text/csv",
            key="download_transcript_screening"
        )

    else:

        st.warning(
            "Transcript screening results are not available."
        )

    st.caption(
        "This analysis uses exact sequence matching against "
        "the retrieved transcript records. It is not genome-wide "
        "specificity analysis."
    )

    # =====================================================
    # STEP 16 — NCBI BLAST SCREENING
    # =====================================================

    st.divider()

    st.subheader(
        "11. NCBI BLAST — Primer Specificity Screening"
    )

    st.warning(
        "This analysis uses the public NCBI BLAST service. "
        "It is not equivalent to Primer-BLAST and does not "
        "provide definitive experimental confirmation."
    )

    organism = st.text_input(
        "Target organism",
        value="Rattus norvegicus",
        help=(
            "Enter the scientific name of the organism."
        ),
        key="rtqpcr_blast_organism"
    )

    render_ncbi_specificity_screen(
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
        organism=organism
    )

    # =====================================================
    # STEP 17 — EXPORT CANDIDATES
    # =====================================================

    st.divider()

    st.subheader("12. Export Primer Candidates")

    candidate_df = pd.DataFrame(candidates)

    export_df = candidate_df.drop(
        columns=["Warnings"],
        errors="ignore"
    )

    csv_data = export_df.to_csv(
        index=False
    ).encode("utf-8-sig")

    st.download_button(
        label="Download All Primer Candidates (CSV)",
        data=csv_data,
        file_name="rtqpcr_primer_candidates.csv",
        mime="text/csv",
        key="download_rtqpcr_candidates"
    )

    # =====================================================
    # STEP 18 — SCIENTIFIC LIMITATIONS
    # =====================================================

    st.divider()

    st.subheader("13. Scientific Interpretation")

    st.markdown(
        """
        **Current analytical limitations:**

        - Primer3 evaluates primer properties and predicted
          primer-pair performance.
        - Hairpin and dimer calculations are computational
          predictions based on thermodynamic models.
        - Transcript screening uses exact sequence matching.
        - Exact transcript matching does not establish
          genome-wide specificity.
        - NCBI BLAST screening is not equivalent to Primer-BLAST.
        - Exon-junction specificity has not yet been implemented.
        - SNP screening has not yet been implemented.
        - Experimental validation is required before using
          any candidate in RT-qPCR.
        """
    )

    st.info(
        "Recommended validation: check promising primer pairs "
        "with NCBI Primer-BLAST, inspect exon structure and "
        "known variants, and experimentally evaluate amplification "
        "efficiency, melt curves, NTC and no-RT controls."
    )