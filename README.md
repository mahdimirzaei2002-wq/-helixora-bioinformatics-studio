# Helixora Bioinformatics Studio

**An academic research workspace for molecular sequence analysis and bioinformatics workflows.**

Helixora brings together gene and transcript retrieval, sequence inspection, coding-variant analysis, public variant database lookup, gene/protein characterization, RT-qPCR primer workflows, and transcript/isoform comparison in one Streamlit application.

## Research modules

1. **Gene & Transcript** — NCBI transcript search, RefSeq retrieval, GenBank annotation and FASTA input.
2. **Sequence Analysis** — sequence inspection, FASTA export and shared transcript registry.
3. **Variant Analysis** — SNV consequence prediction using annotated CDS information.
4. **Database Search** — preliminary ClinVar and dbSNP record retrieval; hits require manual verification.
5. **Transcript & Isoform Analysis** — transcript catalogue, transcript-to-protein mapping, structural comparison and transcript-discriminating region assessment.
6. **Gene & Protein Characterization** — gene/protein annotation, domains, PTMs, expression, literature and structures.
7. **Primer Design** — transcript retrieval, RT-qPCR primer design, candidate evaluation and specificity screening.

## Run locally

Python 3.10–3.12 is recommended.

### Windows PowerShell

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run app.py
```

If PowerShell blocks environment activation:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py
```

### macOS / Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
streamlit run app.py
```

## NCBI configuration

Enter a valid contact email in the **NCBI connection** section of the sidebar. An API key is optional.

For persistent configuration, create `.streamlit/secrets.toml`:

```toml
NCBI_EMAIL = "your.email@example.com"
NCBI_API_KEY = "your_ncbi_api_key"
```

The same settings can be supplied as environment variables. **Never commit API keys, passwords, or private credentials to a public repository.**

For Streamlit Community Cloud, add `NCBI_EMAIL` and (optionally) `NCBI_API_KEY` under the app's Secrets settings. Set the repository's app entry point to `app.py`.

## Scientific notes and limitations

- NCBI, UniProt, GTEx, EBI and BLAST-backed features require internet access and may be affected by external service availability or rate limits.
- Variant annotations are computational predictions and are not clinical classifications.
- ClinVar/dbSNP matches require manual review of the underlying records.
- Transcript variant numbering alone does not establish protein function.
- Transcript-discriminating windows are compared against the retrieved transcript set; this is not genome-wide specificity validation.
- Primer candidates require in-silico specificity assessment and experimental validation.
- Database annotations are reported as supplied by their sources and should not be interpreted as independent experimental confirmation.

**For academic research and educational use only. Not intended for clinical diagnosis or treatment decisions.**

## Project structure

- `app.py` — Streamlit entry point, navigation, session state and presentation layer.
- `config/` — shared application and NCBI configuration.
- `core/` — session-scoped transcript registry.
- `modules/` — scientific analysis and external database integrations.
- `tests/` — automated tests.

## Testing

```bash
python -m compileall -q app.py config modules tests
python -m unittest discover -s tests -v
```

## Developer

**Mahdi Mirzaei**  
Helixora Bioinformatics Studio

## Appearance themes

Helixora includes an in-app light/dark appearance switch in the left sidebar. The selected theme is applied to the app interface, including form fields, placeholders, selectors, code/pattern previews, cards, and data panels. The switch is session-based and resets to light mode when a new session starts.

## Public bilingual website and VPS deployment

The `website/` directory contains the Persian/English project introduction and feature overview. The included Docker + Caddy setup serves the website at `/` and the Streamlit studio at `/app/` on one domain.

See [`DEPLOYMENT.md`](DEPLOYMENT.md) for VPS requirements, DNS, HTTPS and deployment commands. Copy `.env.example` to `.env` and configure your own domain and optional NCBI credentials before deployment.
