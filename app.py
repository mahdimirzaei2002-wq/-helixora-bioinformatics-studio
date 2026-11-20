import streamlit as st
import pandas as pd
import plotly.express as px

from modules.gckr_characterization import render_gckr_page
from modules.primer_design import render_primer_design_page

from modules.ncbi_client import (
    configure_entrez,
    search_transcripts,
    fetch_genbank_record,
    extract_gene_metadata,
    get_cds_feature,
    search_variant_database
)

from modules.mutation_analyzer import analyze_mutation
from modules.mef2c_analysis import render_mef2c_page

from config.settings import (
    AppSettings,
    APP_NAME,
    DEVELOPER_NAME,
    DEVELOPER_EMAIL,
    DEVELOPER_PHONE
)

from core.data_manager import (
    register_active_record,
    get_registry_summary
)


# =====================================================
# 1. PAGE CONFIGURATION
# =====================================================

st.set_page_config(
    page_title="Helixora | Bioinformatics Studio",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


# =====================================================
# 2. GLOBAL STYLE
# =====================================================

# Theme is controlled inside the app so users can switch without changing browser settings.
if "ui_theme" not in st.session_state:
    st.session_state.ui_theme = "light"

DARK = st.session_state.ui_theme == "dark"
# Palette sampled from the supplied Helixora reference code.
# Light mode is a high-contrast counterpart of the same teal/navy identity.
THEME = {
    "canvas": "#0B1220" if DARK else "#F4F7FB",
    "surface": "#14243A" if DARK else "#FFFFFF",
    "surface_alt": "#1B304A" if DARK else "#EAF0F6",
    "ink": "#E8F3FF" if DARK else "#14243A",
    "body": "#C9D8E8" if DARK else "#263B56",
    "muted": "#9FB2C8" if DARK else "#52657A",
    "line": "#263B56" if DARK else "#D5DFEA",
    "input": "#101B2D" if DARK else "#FFFFFF",
    "input_text": "#E8F3FF" if DARK else "#14243A",
    "placeholder": "#9FB2C8" if DARK else "#52657A",
    "brand": "#36D6C5" if DARK else "#087F78",
    "brand_dark": "#20B9AA" if DARK else "#06645F",
    "shadow": "0 12px 34px rgba(0,0,0,.24)" if DARK else "0 12px 34px rgba(20,36,58,.08)",
    "header": "rgba(11,18,32,.94)" if DARK else "rgba(244,247,251,.94)",
    "sidebar": "#101B2D" if DARK else "#EAF0F6",
    "sidebar_ink": "#E8F3FF" if DARK else "#14243A",
    "sidebar_muted": "#9FB2C8" if DARK else "#52657A",
    "sidebar_hover": "rgba(255,255,255,.08)" if DARK else "rgba(8,127,120,.08)",
    "sidebar_selected": "#1B3B50" if DARK else "#D5E9E7",
    "sidebar_selected_border": "#36D6C5" if DARK else "#087F78",
    "card": "#14243A" if DARK else "#FFFFFF",
}

st.markdown(f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Manrope:wght@400;500;600;700;800&display=swap');
:root {{
  --ink:{THEME['ink']}; --body:{THEME['body']}; --muted:{THEME['muted']};
  --canvas:{THEME['canvas']}; --surface:{THEME['surface']}; --surface-alt:{THEME['surface_alt']};
  --line:{THEME['line']}; --brand:{THEME['brand']}; --brand-dark:{THEME['brand_dark']};
  --input:{THEME['input']}; --input-text:{THEME['input_text']}; --placeholder:{THEME['placeholder']};
  --shadow:{THEME['shadow']};
  --sidebar-ink:{THEME['sidebar_ink']}; --sidebar-muted:{THEME['sidebar_muted']};
  --sidebar-hover:{THEME['sidebar_hover']}; --sidebar-selected:{THEME['sidebar_selected']};
  --sidebar-selected-border:{THEME['sidebar_selected_border']};
}}
html, body, [class*="css"] {{ font-family:'DM Sans',sans-serif; }}
.stApp, [data-testid="stAppViewContainer"] {{ background:var(--canvas); color:var(--ink); }}
[data-testid="stHeader"] {{ background:{THEME['header']}; }}
.block-container {{ max-width:1440px; padding-top:1.35rem; padding-bottom:3.2rem; }}
h1,h2,h3,h4 {{ font-family:'Manrope','DM Sans',sans-serif!important; color:var(--ink)!important; letter-spacing:-.035em; }}
h1 {{ font-size:2.05rem!important; font-weight:800!important; }} h2 {{ font-size:1.48rem!important; font-weight:750!important; }} h3 {{ font-size:1.16rem!important; font-weight:750!important; }}
p,li,label,.stMarkdown {{ color:var(--body); }}
[data-testid="stSidebar"] {{ background:{THEME['sidebar']}; border-right:1px solid var(--line); }}
[data-testid="stSidebar"] * {{ color:var(--sidebar-ink); }}
[data-testid="stSidebar"] [data-testid="stRadio"] > label {{ color:var(--sidebar-muted)!important; font-size:.72rem; font-weight:700; letter-spacing:.12em; }}
[data-testid="stSidebar"] [role="radiogroup"] {{ gap:.26rem; }}
[data-testid="stSidebar"] [role="radiogroup"] label {{ padding:.52rem .68rem; border-radius:10px; }}
[data-testid="stSidebar"] [role="radiogroup"] label:hover {{ background:var(--sidebar-hover); }}
[data-testid="stSidebar"] [role="radiogroup"] label:has(input:checked) {{ background:var(--sidebar-selected); border:1px solid var(--sidebar-selected-border); }}
/* Inputs: explicit foreground/background/placeholder colors prevent unreadable email and pattern fields. */
.stTextInput input, .stTextArea textarea, .stNumberInput input,
[data-testid="stSidebar"] [data-testid="stTextInput"] input,
[data-baseweb="input"] input, [data-baseweb="textarea"] textarea {{
  background:var(--input)!important; color:var(--input-text)!important;
  -webkit-text-fill-color:var(--input-text)!important; border:1px solid var(--line)!important;
  border-radius:9px!important; caret-color:var(--brand)!important;
}}
.stTextInput input::placeholder, .stTextArea textarea::placeholder,
[data-baseweb="input"] input::placeholder, [data-baseweb="textarea"] textarea::placeholder,
.stTextInput input::-webkit-input-placeholder, .stTextArea textarea::-webkit-input-placeholder,
.stTextInput input::-moz-placeholder, .stTextArea textarea::-moz-placeholder {{
  color:var(--placeholder)!important; opacity:1!important;
  -webkit-text-fill-color:var(--placeholder)!important; font-weight:600!important;
}}
[data-testid="stSidebar"] [data-testid="stTextInput"] input {{ background:var(--input)!important; color:var(--input-text)!important; -webkit-text-fill-color:var(--input-text)!important; border-color:var(--line)!important; }}
.stTextInput input:focus, .stTextArea textarea:focus, .stNumberInput input:focus,
[data-baseweb="input"] input:focus, [data-baseweb="textarea"] textarea:focus {{ border-color:var(--brand)!important; box-shadow:0 0 0 2px color-mix(in srgb, var(--brand) 22%, transparent)!important; }}
[data-testid="stWidgetLabel"] p, [data-testid="stWidgetLabel"] label {{ color:var(--body)!important; font-weight:600; }}
[data-testid="stSidebar"] [data-testid="stWidgetLabel"] p {{ color:var(--sidebar-ink)!important; }}
[data-testid="stSelectbox"] [data-baseweb="select"] > div,
[data-testid="stMultiSelect"] [data-baseweb="select"] > div {{ background:var(--input)!important; color:var(--input-text)!important; border-color:var(--line)!important; }}
[data-baseweb="popover"], [data-baseweb="menu"], [role="listbox"] {{ background:var(--surface)!important; color:var(--ink)!important; }}
[role="option"] {{ color:var(--ink)!important; }}
[data-testid="stCaptionContainer"] p, .small-muted {{ color:var(--muted)!important; }}
[data-testid="stMetric"] {{ background:var(--surface); border:1px solid var(--line); padding:1rem 1.1rem; border-radius:15px; box-shadow:var(--shadow); }}
[data-testid="stMetricLabel"] {{ color:var(--muted)!important; font-weight:600; }}
[data-testid="stMetricValue"] {{ color:var(--ink)!important; font-family:'Manrope',sans-serif; font-weight:800; }}
div.stButton>button,div.stDownloadButton>button {{ border-radius:10px; min-height:2.65rem; font-weight:700; border:1px solid var(--line); background:var(--surface); color:var(--ink); transition:all .18s ease; }}
div.stButton>button[kind="primary"] {{ background:var(--brand); border-color:var(--brand); color:white; }}
div.stButton>button[kind="primary"]:hover {{ background:var(--brand-dark); border-color:var(--brand-dark); transform:translateY(-1px); }}
div.stButton>button:hover,div.stDownloadButton>button:hover {{ border-color:var(--brand); color:var(--brand); }}
div[data-testid="stExpander"] {{ background:var(--surface); border:1px solid var(--line); border-radius:13px; }}
[data-testid="stExpander"] summary, [data-testid="stExpander"] summary p {{ color:var(--ink)!important; }}
[data-testid="stDataFrame"], [data-testid="stTable"] {{ border:1px solid var(--line); border-radius:12px; overflow:hidden; }}
.stAlert {{ border-radius:12px; }}
hr {{ border-color:var(--line); }}
.brand-lockup {{ display:flex; align-items:center; gap:12px; padding:.25rem 0 .9rem; }}
.brand-mark {{ width:42px;height:42px;border-radius:13px;display:grid;place-items:center;background:linear-gradient(135deg,#2AA895,#126B68);color:white;font-size:22px;box-shadow:0 7px 18px rgba(42,168,149,.22); }}
.brand-name {{ font:800 1.16rem 'Manrope',sans-serif; letter-spacing:-.04em; color:var(--sidebar-ink); }}
.brand-tag {{ color:var(--sidebar-muted); font-size:.68rem; letter-spacing:.12em; text-transform:uppercase; margin-top:2px; }}
.topline {{ display:flex;justify-content:space-between;align-items:center;padding:.15rem 0 1rem;border-bottom:1px solid var(--line);margin-bottom:1.5rem; }}
.topline-brand {{ color:var(--brand);font:800 .92rem 'Manrope',sans-serif;letter-spacing:-.02em; }}
.topline-label {{ color:var(--muted);font-size:.76rem;font-weight:700;letter-spacing:.1em;text-transform:uppercase; }}
.hero-panel {{ padding:1.7rem 1.8rem;border-radius:20px;background:linear-gradient(115deg,#102F40 0%,#155B60 62%,#218778 100%);color:white;box-shadow:0 16px 36px rgba(16,47,64,.15);margin:.35rem 0 1.25rem; }}
.hero-eyebrow {{ color:#A8E6D9;text-transform:uppercase;letter-spacing:.16em;font-size:.7rem;font-weight:800; }}
.hero-title {{ color:white;font:800 2rem 'Manrope',sans-serif;letter-spacing:-.05em;margin:.45rem 0; }}
.hero-copy {{ color:#E0F0F0;font-size:.96rem;max-width:760px;line-height:1.65; }}
.module-card {{ background:var(--surface);border:1px solid var(--line);border-radius:15px;padding:1rem 1.05rem;min-height:128px;box-shadow:var(--shadow); }}
.module-card:hover {{ border-color:var(--brand); }}
.module-icon {{ font-size:1.25rem;margin-bottom:.35rem; }}
.module-title {{ color:var(--ink);font:750 .98rem 'Manrope',sans-serif;margin-bottom:.28rem; }}
.module-copy {{ color:var(--muted);font-size:.82rem;line-height:1.5; }}
.section-kicker {{ color:var(--brand);font-size:.7rem;font-weight:800;letter-spacing:.14em;text-transform:uppercase;margin-bottom:.3rem; }}
/* Code, sequence previews and pattern-like panels */
code, pre, [data-testid="stCode"] {{ color:{'#DCEAF5' if DARK else '#17324A'}!important; background:{'#102033' if DARK else '#EAF1F7'}!important; border:1px solid var(--line); border-radius:8px; }}
[data-testid="stMarkdownContainer"] table {{ color:var(--ink); }}
@media (max-width: 700px) {{ .block-container {{ padding-left:1rem;padding-right:1rem; }} .hero-title {{ font-size:1.55rem; }} .topline-label {{ display:none; }} }}
</style>
""", unsafe_allow_html=True)


# =====================================================
# 3. SESSION STATE
# =====================================================

def init_state():

    defaults = {

        "record": None,
        "cds_info": None,
        "metadata": None,

        "analysis_df": None,
        "db_results": None,

        "selected_accession": None,

        "gene_candidates": [],
        "gene_symbol": "",
        "organism": "Homo sapiens",

        "ncbi_email": "",
        "ncbi_api_key": "",

        "primer_transcripts": [],

        "transcript_isoform_records": []

    }

    for key, value in defaults.items():

        if key not in st.session_state:

            st.session_state[key] = value


init_state()


# =====================================================
# 4. NCBI CENTRAL CONFIGURATION
# =====================================================

central_settings = AppSettings.load()

if not st.session_state.ncbi_email:

    st.session_state.ncbi_email = central_settings.ncbi_email

if not st.session_state.ncbi_api_key:

    st.session_state.ncbi_api_key = central_settings.ncbi_api_key


# =====================================================
# 5. NCBI CONNECTION
# =====================================================

def ncbi_ready():

    email = st.session_state.ncbi_email.strip()

    api_key = st.session_state.ncbi_api_key.strip()

    try:

        configure_entrez(
            email=email,
            api_key=api_key
        )

        return True

    except Exception as exc:

        st.error(
            f"NCBI configuration error: {exc}"
        )

        return False


# =====================================================
# 6. SIDEBAR NAVIGATION
# =====================================================

with st.sidebar:
    st.markdown("""
    <div class="brand-lockup">
      <div class="brand-mark">⌬</div>
      <div><div class="brand-name">HELIXORA</div><div class="brand-tag">Bioinformatics Studio</div></div>
    </div>
    """, unsafe_allow_html=True)
    st.markdown(f'<div style="color:{THEME['sidebar_muted']};font-size:.78rem;line-height:1.5">An academic workspace for sequence intelligence and molecular research.</div>', unsafe_allow_html=True)
    theme_label = "☀️  Switch to light mode" if DARK else "🌙  Switch to dark mode"
    if st.button(theme_label, key="theme_toggle", use_container_width=True):
        st.session_state.ui_theme = "light" if DARK else "dark"
        st.rerun()
    st.divider()
    page = st.radio(
        "RESEARCH WORKSPACE",
        [
            "Overview",
            "Gene & Transcript",
            "Sequence Analysis",
            "Variant Analysis",
            "Database Search",
            "Transcript & Isoform Analysis",
            "Gene & Protein Characterization",
            "Primer Design",
            "Settings & About"
        ],
        key="main_navigation"
    )
    st.divider()
    with st.expander("NCBI connection", expanded=False):
        st.caption("Shared credentials for NCBI-enabled workflows.")
        st.text_input("Contact email", key="ncbi_email", help="Used for shared Entrez requests.")
        st.text_input("API key (optional)", type="password", key="ncbi_api_key", help="Optional; increases the permitted request rate.")
        st.caption("An internet connection is required for live database retrieval.")
    st.markdown('<div style="padding-top:1.2rem;color:#8EABB5;font-size:.72rem;line-height:1.6">HELIXORA · ACADEMIC EDITION<br>Research tools, not clinical diagnostics.</div>', unsafe_allow_html=True)


# =====================================================
# 7. MAIN HEADER
# =====================================================

st.markdown(f"""
<div class="topline">
  <div class="topline-brand">HELIXORA <span style="color:#91A6B5;font-weight:500">/ Molecular Research Workspace</span></div>
  <div class="topline-label">{page}</div>
</div>
""", unsafe_allow_html=True)


# =====================================================
# PAGE 1 — OVERVIEW
# =====================================================

if page == "Overview":
    record = st.session_state.record
    analysis_df = st.session_state.analysis_df
    st.markdown("""
    <div class="hero-panel">
      <div class="hero-eyebrow">Academic bioinformatics environment</div>
      <div class="hero-title">Explore biology. Decode the sequence.</div>
      <div class="hero-copy">A connected research workspace for gene discovery, transcript annotation, sequence analysis, variant interpretation, protein characterization and molecular assay design.</div>
    </div>
    """, unsafe_allow_html=True)

    c1, c2, c3 = st.columns(3)
    c1.metric("Active transcript", record.id if record else "Not loaded")
    c2.metric("Sequence length", f"{len(record.seq):,} nt" if record else "—")
    c3.metric("Variants analyzed", len(analysis_df) if analysis_df is not None else 0)

    st.markdown('<div class="section-kicker">Your research toolkit</div>', unsafe_allow_html=True)
    st.markdown("## Integrated analysis modules")
    modules = [
        ("🧬", "Gene & Transcript", "Retrieve NCBI records, explore RefSeq accessions and inspect coding annotations."),
        ("⌁", "Sequence Analysis", "Inspect nucleotide sequences, export FASTA and access the shared transcript registry."),
        ("◎", "Variant Analysis", "Evaluate substitutions and explore predicted coding consequences and HGVS notation."),
        ("◈", "Database Search", "Retrieve preliminary ClinVar and dbSNP records for further review."),
        ("⌘", "Transcript & Isoform Analysis", "Compare transcript structures, map proteins and examine alternative-splicing evidence."),
        ("◉", "Gene & Protein Characterization", "Explore gene annotations, protein domains, PTMs, expression and structures."),
        ("✳", "Primer Design", "Develop and assess molecular assays, including RT-qPCR primer workflows.")
    ]
    for i in range(0, len(modules), 2):
        left, right = st.columns(2)
        for col, item in zip((left, right), modules[i:i+2]):
            icon, title, description = item
            with col:
                st.markdown(f'<div class="module-card"><div class="module-icon">{icon}</div><div class="module-title">{title}</div><div class="module-copy">{description}</div></div>', unsafe_allow_html=True)
    st.markdown("")
    with st.expander("Suggested research workflow", expanded=False):
        st.markdown("""
        1. Identify a target gene and organism.
        2. Retrieve transcript records and review annotations.
        3. Inspect sequence and coding-region information.
        4. Analyze variants and consult public databases.
        5. Compare transcript and protein isoforms.
        6. Design and evaluate molecular assays.
        7. Export results for downstream research.
        """)
    st.caption("Scientific outputs are computational aids. Review source annotations and validate findings before drawing biological conclusions.")


# =====================================================
# PAGE 2 — GENE & TRANSCRIPT
# =====================================================

elif page == "Gene & Transcript":

    st.subheader("Find and load a transcript")

    mode = st.radio(

        "Sequence identification method",

        [

            "Gene symbol",

            "RefSeq accession",

            "Paste FASTA"

        ],

        horizontal=True

    )


    # -----------------------------------------------
    # GENE SYMBOL SEARCH
    # -----------------------------------------------

    if mode == "Gene symbol":

        col1, col2 = st.columns([2, 1])

        with col1:

            gene_symbol = st.text_input(

                "Gene symbol",

                value=st.session_state.gene_symbol,

                placeholder="e.g. TP53"

            )

        with col2:

            organism = st.text_input(

                "Organism",

                value=st.session_state.organism

            )

        if st.button(

            "Search NCBI Transcripts",

            type="primary"

        ):

            if gene_symbol.strip() and ncbi_ready():

                try:

                    with st.spinner(
                        "Searching NCBI nuccore..."
                    ):

                        candidates = search_transcripts(

                            gene_symbol.strip(),

                            organism.strip()

                        )

                    st.session_state.gene_candidates = candidates

                    st.session_state.gene_symbol = gene_symbol.strip()

                    st.session_state.organism = organism.strip()

                    if candidates:

                        st.success(
                            f"{len(candidates)} transcript candidates found."
                        )

                    else:

                        st.warning(
                            "No transcript candidates were found."
                        )

                except Exception as exc:

                    st.error(
                        f"NCBI search failed: {exc}"
                    )

        if st.session_state.gene_candidates:

            labels = [

                item["label"]

                for item in st.session_state.gene_candidates

            ]

            chosen = st.selectbox(

                "Select transcript",

                labels

            )

            selected = st.session_state.gene_candidates[
                labels.index(chosen)
            ]

            st.caption(
                "Review transcript accession and description before loading."
            )

            if st.button(

                "Download GenBank Record",

                type="primary"

            ):

                if ncbi_ready():

                    try:

                        with st.spinner(
                            f"Downloading {selected['accession']}..."
                        ):

                            record = fetch_genbank_record(
                                selected["accession"]
                            )

                            cds_info = get_cds_feature(record)

                            metadata = extract_gene_metadata(record)

                        st.session_state.record = record

                        st.session_state.cds_info = cds_info

                        st.session_state.metadata = metadata

                        st.session_state.selected_accession = record.id

                        st.session_state.analysis_df = None

                        st.session_state.db_results = None

                        register_active_record(

                            record,

                            gene=st.session_state.gene_symbol,

                            organism=st.session_state.organism

                        )

                        st.success(
                            "Transcript and annotation loaded successfully."
                        )

                    except Exception as exc:

                        st.error(
                            f"Record retrieval failed: {exc}"
                        )


    # -----------------------------------------------
    # REFSEQ ACCESSION
    # -----------------------------------------------

    elif mode == "RefSeq accession":

        accession = st.text_input(

            "RefSeq accession",

            placeholder="e.g. NM_003106.4"

        )

        if st.button(

            "Download Transcript",

            type="primary"

        ):

            if accession.strip() and ncbi_ready():

                try:

                    with st.spinner(
                        "Retrieving GenBank record..."
                    ):

                        record = fetch_genbank_record(
                            accession.strip()
                        )

                        cds_info = get_cds_feature(record)

                        metadata = extract_gene_metadata(record)

                    st.session_state.record = record

                    st.session_state.cds_info = cds_info

                    st.session_state.metadata = metadata

                    st.session_state.selected_accession = record.id

                    st.session_state.analysis_df = None

                    st.session_state.db_results = None

                    register_active_record(record)

                    st.success(
                        "Transcript downloaded successfully."
                    )

                except Exception as exc:

                    st.error(
                        f"Transcript retrieval failed: {exc}"
                    )


    # -----------------------------------------------
    # FASTA INPUT
    # -----------------------------------------------

    else:

        uploaded = st.file_uploader(

            "Upload FASTA",

            type=["fa", "fasta", "fna", "txt"]

        )

        pasted = st.text_area(

            "Or paste FASTA sequence",

            height=180,

            placeholder=">Transcript_ID\nATG..."

        )

        if st.button(

            "Load FASTA",

            type="primary"

        ):

            try:

                from Bio import SeqIO

                from io import StringIO

                if pasted.strip():

                    record = SeqIO.read(
                        StringIO(pasted),
                        "fasta"
                    )

                elif uploaded:

                    record = SeqIO.read(

                        StringIO(
                            uploaded.getvalue().decode("utf-8")
                        ),

                        "fasta"

                    )

                else:

                    raise ValueError(
                        "Upload or paste a FASTA sequence."
                    )

                st.session_state.record = record

                st.session_state.metadata = extract_gene_metadata(record)

                st.session_state.cds_info = None

                st.session_state.selected_accession = record.id

                st.session_state.analysis_df = None

                st.session_state.db_results = None

                register_active_record(record)

                st.success("FASTA sequence loaded.")

            except Exception as exc:

                st.error(
                    f"FASTA loading failed: {exc}"
                )


    # -----------------------------------------------
    # TRANSCRIPT INFORMATION
    # -----------------------------------------------

    if st.session_state.record is not None:

        st.divider()

        st.markdown("### Transcript record")

        record = st.session_state.record

        st.success(
            f"Active transcript: {record.id}"
        )

        left, right = st.columns(2)

        with left:

            st.markdown("#### Record Metadata")

            st.json(
                st.session_state.metadata or {}
            )

        with right:

            st.markdown("#### CDS Annotation")

            if st.session_state.cds_info:

                cds = st.session_state.cds_info

                st.metric(
                    "CDS Start",
                    cds["start"]
                )

                st.metric(
                    "CDS End",
                    cds["end"]
                )

                st.metric(
                    "CDS Length",
                    f"{len(cds['sequence']):,} nt"
                )

            else:

                st.warning(
                    "No CDS annotation was detected."
                )

                c1, c2 = st.columns(2)

                manual_start = c1.number_input(

                    "Manual CDS start",

                    min_value=1,

                    value=1,

                    key="manual_cds_start"

                )

                manual_end = c2.number_input(

                    "Manual CDS end",

                    min_value=1,

                    value=len(record.seq),

                    key="manual_cds_end"

                )

                if st.button("Confirm CDS Boundaries"):

                    if manual_start <= manual_end <= len(record.seq):

                        st.session_state.cds_info = {

                            "start": int(manual_start),

                            "end": int(manual_end),

                            "sequence": str(
                                record.seq[
                                    int(manual_start)-1:int(manual_end)
                                ]
                            ),

                            "feature": None

                        }

                        st.success(
                            "Manual CDS boundaries saved."
                        )

                    else:

                        st.error(
                            "Invalid CDS boundaries."
                        )


# =====================================================
# PAGE 3 — SEQUENCE ANALYSIS
# =====================================================

elif page == "Sequence Analysis":

    st.subheader("Sequence Analysis")

    record = st.session_state.record

    if record is None:

        st.info(
            "Load a transcript in Gene & Transcript first."
        )

    else:

        st.metric(
            "Sequence Length",
            f"{len(record.seq):,} nt"
        )

        st.code(
            str(record.seq),
            language="text"
        )

        st.download_button(

            "Download FASTA",

            f">{record.id}\n{str(record.seq)}\n",

            file_name=f"{record.id}.fasta",

            mime="text/plain"

        )

        st.markdown("#### Shared Transcript Registry")

        shared = get_registry_summary()

        if shared:

            st.dataframe(

                pd.DataFrame(shared),

                use_container_width=True,

                hide_index=True

            )


# =====================================================
# PAGE 4 — VARIANT ANALYSIS
# =====================================================

elif page == "Variant Analysis":

    st.subheader("Variant analysis")

    record = st.session_state.record

    cds_info = st.session_state.cds_info

    if record is None:

        st.warning(
            "Load a transcript first."
        )

    elif cds_info is None:

        st.warning(
            "A CDS annotation is required."
        )

    else:

        st.info(

            f"Analyzing {record.id} | "

            f"CDS: {cds_info['start']}–{cds_info['end']}"

        )

        coordinate_system = st.radio(

            "Coordinate system",

            ["mRNA", "CDS"],

            horizontal=True

        )

        n = st.number_input(

            "Number of variants",

            min_value=1,

            max_value=100,

            value=1

        )

        variant_inputs = []

        for i in range(int(n)):

            with st.expander(
                f"Variant {i+1}",
                expanded=(i == 0)
            ):

                a, b, c = st.columns(3)

                pos = a.number_input(

                    "Position",

                    min_value=1,

                    value=1,

                    key=f"vpos_{i}"

                )

                ref = b.selectbox(

                    "Reference",

                    ["A", "C", "G", "T"],

                    key=f"vref_{i}"

                )

                alt = c.selectbox(

                    "Alternate",

                    ["A", "C", "G", "T"],

                    key=f"valt_{i}"

                )

                variant_inputs.append(
                    (int(pos), ref, alt)
                )

        if st.button(

            "Run Variant Analysis",

            type="primary",

            use_container_width=True

        ):

            rows = []

            errors = []

            for i, (pos, ref, alt) in enumerate(

                variant_inputs,

                start=1

            ):

                try:

                    row = analyze_mutation(

                        sequence=record.seq,

                        cds_start=cds_info["start"],

                        cds_end=cds_info["end"],

                        position=pos,

                        reference=ref,

                        alternate=alt,

                        coordinate_system=coordinate_system,

                        accession=record.id

                    )

                    row["Variant #"] = i

                    rows.append(row)

                except Exception as exc:

                    errors.append(
                        f"Variant {i}: {exc}"
                    )

            for error in errors:

                st.error(error)

            if rows:

                st.session_state.analysis_df = pd.DataFrame(rows)

                st.session_state.db_results = None

                st.success(
                    f"{len(rows)} variant(s) analyzed."
                )

        if st.session_state.analysis_df is not None:

            df = st.session_state.analysis_df

            st.markdown("### Results")

            st.dataframe(

                df,

                use_container_width=True,

                hide_index=True

            )

            if "Consequence" in df.columns:

                metrics = st.columns(4)

                categories = [

                    "Synonymous",

                    "Missense",

                    "Stop gained",

                    "Stop lost"

                ]

                for col, category in zip(metrics, categories):

                    col.metric(

                        category,

                        int(
                            (df["Consequence"] == category).sum()
                        )

                    )

                chart = px.pie(

                    df,

                    names="Consequence",

                    title="Predicted Coding Consequences",

                    hole=0.55

                )

                st.plotly_chart(

                    chart,

                    use_container_width=True

                )

            st.download_button(

                "Download Variant Report",

                df.to_csv(
                    index=False
                ).encode("utf-8-sig"),

                file_name="variant_analysis.csv",

                mime="text/csv"

            )


# =====================================================
# PAGE 5 — DATABASE SEARCH
# =====================================================

elif page == "Database Search":

    st.subheader("Public variant database search")

    df = st.session_state.analysis_df

    if df is None or df.empty:

        st.warning(
            "Run variant analysis first."
        )

    else:

        st.warning(
            "Database search results require manual verification."
        )

        st.dataframe(

            df,

            use_container_width=True,

            hide_index=True

        )

        if st.button(

            "Search ClinVar and dbSNP",

            type="primary"

        ):

            if ncbi_ready():

                combined = {}

                progress = st.progress(0)

                hgvs_values = df[
                    "HGVS (transcript)"
                ].tolist()

                for i, hgvs in enumerate(hgvs_values):

                    try:

                        combined[hgvs] = search_variant_database(hgvs)

                    except Exception as exc:

                        combined[hgvs] = {
                            "error": str(exc)
                        }

                    progress.progress(
                        (i + 1) / len(hgvs_values)
                    )

                st.session_state.db_results = combined

                st.success(
                    "Database searches completed."
                )

        if st.session_state.db_results:

            for hgvs, result in st.session_state.db_results.items():

                with st.expander(hgvs):

                    if "error" in result:

                        st.error(result["error"])

                        continue

                    for db in ("clinvar", "snp"):

                        item = result.get(db, {})

                        st.markdown(
                            f"**{db.upper()}**"
                        )

                        st.write(
                            "Status:",
                            item.get("status", "No result")
                        )

                        st.write(
                            "Record count:",
                            item.get("count", 0)
                        )

                        if item.get("ids"):

                            st.write(
                                "NCBI IDs:",
                                ", ".join(item["ids"])
                            )

                        if item.get("summaries"):

                            st.dataframe(

                                pd.DataFrame(item["summaries"]),

                                use_container_width=True,

                                hide_index=True

                            )


# =====================================================
# PAGE 6 — GENE & PROTEIN CHARACTERIZATION
# =====================================================

elif page == "Gene & Protein Characterization":

    render_gckr_page()


# =====================================================
# PAGE 7 — PRIMER DESIGN
# =====================================================

elif page == "Primer Design":

    if ncbi_ready():

        render_primer_design_page()


# =====================================================
# PAGE 8 — TRANSCRIPT & ISOFORM ANALYSIS
# =====================================================

elif page == "Transcript & Isoform Analysis":

    if ncbi_ready():

        render_mef2c_page()


# =====================================================
# PAGE 9 — SETTINGS & ABOUT
# =====================================================

elif page == "Settings & About":

    st.subheader(
        "Application Settings and Project Information"
    )

    st.markdown("### Shared NCBI Configuration")

    st.write(
        "Configure the NCBI contact email and optional API key "
        "in the sidebar. All NCBI-enabled modules use the "
        "shared Entrez configuration."
    )

    st.code(
        "NCBI_EMAIL=your.email@example.com\n"
        "NCBI_API_KEY=your_api_key",
        language="ini"
    )

    st.caption(
        "For persistent configuration, use Streamlit secrets "
        "or environment variables. Never commit API keys "
        "to source control."
    )

    st.divider()

    st.markdown("### About This Application")

    st.markdown(f"""
    **Helixora Bioinformatics Studio**

    An integrated academic research platform for molecular
    biology and bioinformatics workflows.

    **Developer:** Mahdi Mirzaei
### Research Modules

    **1. Gene & Transcript**
    - NCBI transcript retrieval
    - RefSeq accession search
    - GenBank annotation
    - CDS identification
    - FASTA input

    **2. Sequence Analysis**
    - Sequence inspection
    - FASTA export
    - Shared transcript registry

    **3. Variant Analysis**
    - SNV analysis
    - Coding consequence prediction
    - HGVS generation

    **4. Database Search**
    - ClinVar retrieval
    - dbSNP retrieval

    **5. Gene & Protein Characterization**
    - Gene and protein annotation
    - SNP information
    - Protein domains
    - Post-translational modifications
    - Literature
    - Protein structures

    **6. Primer Design**
    - Transcript retrieval
    - RT-qPCR primer design
    - Primer evaluation
    - Specificity screening

    **7. Transcript & Isoform Analysis**
    - Gene-independent transcript analysis
    - Transcript-to-protein mapping
    - Alternative splicing analysis
    - Transcript comparison
    - Isoform-specific primer feasibility

    ### Important Scientific Limitations

    - Transcript variant numbering does not establish protein function.
    - Functional isoform classification requires supporting evidence.
    - Predicted mutation consequences are not clinical classifications.
    - Primer candidates require specificity assessment.
    - Experimental validation is necessary before biological interpretation.

    This platform is intended for academic research and educational use,
    not clinical diagnosis.
    """)

    st.divider()

    st.caption(
        "Helixora Bioinformatics Studio | Academic Research Platform"
    )
